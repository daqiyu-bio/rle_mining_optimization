#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------------------#
# Copyright (c) 2026 Beijing Astragenomics Technology Co., Ltd. All rights reserved.        #
# Licensing: see the LICENSE file in the repository root.                                   #
# version : 1.0                                                                             #
#                                                                                           #
# Author:                                                                                   #
#   Kai Liu <kai.liu@astragenomics.com>                                                     #
#                                                                                           #
# Description:                                                                              #
#   Annotate insertion sites with gene-relative position, chromatin state, and              #
#   safe harbor features.                                                                   #
# ------------------------------------------------------------------------------------------#

import os
import tqdm
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from multiprocessing import Pool
from subprocess import getstatusoutput
import time


def parse_params():
    parser = ArgumentParser()
    parser.add_argument('--mapping', dest='mapping_file', action='store', help='mapping result site, use 0-base coordinate.', required=True)
    parser.add_argument('--gtf', dest='gtf_file', action='store', help='GTF file from UCSC trimmed by table format.', required=True)  # GTF gene body coordinates are 1-based in this pipeline
    parser.add_argument('--cancergene', dest='cancer_genelist', action='store', help='Cancer gene list from OncoKB.', required=True)
    parser.add_argument('--conserved', dest='conserved_region', action='store', help='conserved region in genome from UCSC.', required=True)
    parser.add_argument('--miRNA', dest='miRNA', action='store', help='miRNA gene region from UCSC.', required=True)
    parser.add_argument('--chrom_state', dest='chrom_state', action='store', help='chromatin states from epigenomics project.', required=True)
    parser.add_argument('--outpath', dest='outpath', action='store', help='specify the output directory.', required=True)
    parser.add_argument('--outlabel', dest='outlabel', action='store', help='specify the output file label.', required=True)
    parser.add_argument('--thread', dest='thread', action='store', help='Number of worker processes for multiprocessing pools.', default=10, type=int)
    parser.add_argument('--fai',dest='fai',action='store',help='related genome.fai file(samtools faidx..)',type=str)
    params = parser.parse_args()
    return params


def cal_distance(start, end, start_, end_, strand):
    """
    Calculate the signed distance and relative orientation between an insertion site
    and a gene body.
    """
    center_distance = np.mean([start, end]) - np.mean([start_, end_])
    if (center_distance >= 0 and strand == '+') or (center_distance < 0 and strand == '-'):
        relative_pos = 'downstream'  
    else:
        relative_pos = 'upstream'    

    distance = np.abs(center_distance) - (end - start) / 2 - (end_ - start_) / 2
    if distance < 0: 
        relative_pos = 'inner'
    else:
        pass
    return distance, relative_pos    

def get_TSS(x):
    if x['strand'] == '-':
        return x['end'], x['end'] + 1
    else:
        return x['start'] - 1, x['start']


def cal_distance_groupby_chromosome(insert_site, refer_info_, step=5000):
    """
    For each insertion site, expand a genomic window until genes are found and pick the nearest hit.
    """
    if len(refer_info_) == 0:
        return insert_site

    for i in tqdm.tqdm(insert_site.index):
        start = insert_site.loc[i, 'start']
        end = insert_site.loc[i, 'end']

        expand = 1
        while refer_info_.query('(start > @start - @step * @expand and end < @end + @step * @expand) or \
                                (start < @start - @step * @expand  and end > @start - @step * @expand) or \
                                (start < @end + @step * @expand and end > @end + @step * @expand)').shape[0] < 1:
            expand *= 2
        refer_info = refer_info_.query('(start > @start - @step * @expand and end < @end + @step * @expand) or \
                                (start < @start - @step * @expand  and end > @start - @step * @expand) or \
                                (start < @end + @step * @expand and end > @end + @step * @expand)').copy()

        for j in refer_info.index:
            start_ = refer_info.loc[j, 'start']
            end_ = refer_info.loc[j, 'end']
            strand = refer_info.loc[j, 'strand']
            distance, relative_pos = cal_distance(start, end, start_, end_, strand)
            refer_info.loc[j, 'distance'] = distance
            refer_info.loc[j, 'relative_pos_gene'] = relative_pos
        sele = ['gene_name', 'distance', 'relative_pos_gene']
        insert_site.loc[i, sele] = refer_info.sort_values(['distance']).head(1)[sele].iloc[0]
        if insert_site.loc[i,'relative_pos_gene'] == 'inner':
            inner_genes = refer_info.query('distance<0').drop_duplicates('gene_name')['gene_name'].to_list()
            insert_site.loc[i, 'gene_name'] = ';'.join(inner_genes)
    return insert_site


def grant_class_by_gene_distance(x):
    if x['relative_pos_gene'] == 'inner':
        return '0.0kb'
    else:
        level = x['distance'] // 10000
        if level <= 4:
            return f'{level+1:.0f}0kb'
        else:
            return '>50kb'


def grant_class_by_TSS_distance(x):
    """
    Classify distance to TSS: distances ≤250 nt are assigned to the TSS category(0.0kb);
    larger distances are binned into 500 nt intervals up to >2.5 kb.
    """
    if x['distance'] <= 250:
        return '0.0kb'
    else:
        pos = (x['distance'] - 250) // 500
        if pos <= 4:
            return f'{0.5 * (pos +1)}kb'
        else:
            return '>2.5kb'

def grant_class_by_cancerGene(x):
    if x['distance'] <= 1000 and x['relative_pos_gene'] in ['upstream','inner']:
        return 'Yes'
    else:
        return 'No'


def get_step(chrom,candidate_num):
    if candidate_num <= 0:
        candidate_num = 1
    if chrom in CHROM_LENGTH:
        step = CHROM_LENGTH[chrom]/candidate_num
    else:
        step = 50000
    return int(step)


if __name__ == '__main__':
    start_time = time.time()
    params = parse_params()
    mapping_file = params.mapping_file
    gtf_file = params.gtf_file
    cancer_genelist = params.cancer_genelist
    conserved_region_file = params.conserved_region
    miRNA_file = params.miRNA
    chrom_state = params.chrom_state
    outpath = params.outpath
    outlabel = params.outlabel
    thread = params.thread

    CHROM_LENGTH = pd.read_csv(params.fai ,sep='\t',header=None).set_index(0)[1].to_dict()
    os.makedirs(os.path.join(outpath, 'temp'), exist_ok=True)
    os.makedirs(os.path.join(outpath, 'stat'), exist_ok=True)

    mapping_site = pd.read_csv(mapping_file, sep='\t', header=None)
    if len(mapping_site.columns) > 4:
        print('Mapping input must be 3–4 column BED without a header; extra columns are ignored.')
        mapping_site = mapping_site.iloc[:, 0:4].copy()  
    mapping_site.columns = ['chrom', 'start', 'end', 'other'][:min(4, len(mapping_site.columns))]  
    mapping_site['pos_label'] = mapping_site.apply(lambda x: f'{x["chrom"]}_{x["start"]}_{x["end"]}', axis=1)
    mapping_site_ = mapping_site.copy()

    all_gene_Info = pd.read_csv(gtf_file, sep='\t',low_memory=False) 
    all_gene_Info_transcript = all_gene_Info.query('Level=="transcript"').copy()
    all_gene_Info_TSS = all_gene_Info_transcript.copy()

    all_gene_Info_transcript['start'] = all_gene_Info_transcript['start'] - 1 
    all_gene_Info_TSS[['start', 'end']] = all_gene_Info_TSS.apply(get_TSS, axis=1, result_type='expand')  

    work_space = []
    pool = Pool(thread)
    for chrom in mapping_site_['chrom'].drop_duplicates():
        print(chrom)
        insert_site = mapping_site_.query('chrom == @chrom').copy()
        refer_info_ = all_gene_Info_transcript.query('chrom == @chrom')
        insert_sites = np.array_split(insert_site, np.ceil(len(insert_site) / 5000))
        work_space.extend([pool.apply_async(cal_distance_groupby_chromosome, args=(i, refer_info_, 50000)) for i in insert_sites])
    pool.close()
    pool.join()
    mapping_site_gene_distance = pd.concat([i.get() for i in work_space])
    mapping_site_gene_distance['gene_distance_class'] = mapping_site_gene_distance.apply(grant_class_by_gene_distance, axis=1)

    mapping_site_gene_distance['out_gene'] = mapping_site_gene_distance.apply(lambda x: 'No' if x['relative_pos_gene'] == 'inner' else 'Yes', axis=1)

    all_gene_Info_Exon = all_gene_Info.query('Level=="exon"').copy()
    work_space = []
    pool = Pool(thread)
    for chrom in mapping_site_['chrom'].drop_duplicates():
        insert_site = mapping_site_.query('chrom == @chrom').copy()
        refer_info_ = all_gene_Info_Exon.query('chrom == @chrom')
        insert_sites = np.array_split(insert_site, np.ceil(len(insert_site) / 5000))
        work_space.extend([pool.apply_async(cal_distance_groupby_chromosome, args=(i, refer_info_, get_step(chrom,len(refer_info_)))) for i in insert_sites])
    pool.close()
    pool.join()
    mapping_site_Exon_distance = pd.concat([i.get() for i in work_space])
    mapping_site_Exon_distance.to_csv(os.path.join(outpath, 'temp', f'{outlabel}.relative_pos_Exon.csv'), sep='\t', index=False)
    In_Exon_pos = mapping_site_Exon_distance.query('relative_pos_gene=="inner"')['pos_label'].to_list()

    all_gene_Info_Intron = all_gene_Info.query('Level=="intron"').copy()
    work_space = []
    pool = Pool(thread)
    for chrom in mapping_site_['chrom'].drop_duplicates():
        insert_site = mapping_site_.query('chrom == @chrom and ~pos_label.isin(@In_Exon_pos)').copy()
        if len(insert_site) == 0:
            continue
        refer_info_ = all_gene_Info_Intron.query('chrom == @chrom')
        try:
            insert_sites = np.array_split(insert_site, np.ceil(len(insert_site) / 5000))
        except Exception as error:
            print('Failed to split tasks...')
            raise RuntimeError('Failed to split insert_site chunks') from error
        work_space.extend([pool.apply_async(cal_distance_groupby_chromosome, args=(i, refer_info_, get_step(chrom,len(refer_info_)))) for i in insert_sites])
    pool.close()
    pool.join()
    mapping_site_Intron_distance = pd.concat([i.get() for i in work_space])
    mapping_site_Intron_distance.to_csv(os.path.join(outpath, 'temp', f'{outlabel}.relative_pos_Intron.csv'), sep='\t', index=False)
    In_Intron_pos = mapping_site_Intron_distance.query('relative_pos_gene=="inner"')['pos_label'].to_list()

    mapping_site_gene_distance.loc[mapping_site_gene_distance['pos_label'].isin(In_Exon_pos), 'In_gene_type'] = 'Exon'
    mapping_site_gene_distance.loc[mapping_site_gene_distance['pos_label'].isin(In_Intron_pos), 'In_gene_type'] = 'Intron'
    mapping_site_gene_distance['In_gene_type'] = mapping_site_gene_distance['In_gene_type'].fillna('out_gene')

    mapping_site = pd.merge(mapping_site, mapping_site_gene_distance[['pos_label', 'relative_pos_gene', 'gene_distance_class', 'out_gene', 'In_gene_type']].drop_duplicates('pos_label'), left_on=['pos_label'], right_on=['pos_label'], how='left')

    mapping_site_gene_distance.to_csv(os.path.join(outpath, 'temp', f'{outlabel}.relative_pos_gene.csv'), sep='\t', index=False)
    mapping_site.to_csv(os.path.join(outpath, f'{outlabel}.1.csv'), sep='\t', index=False)

    work_space = []
    pool = Pool(thread)
    for chrom in mapping_site_['chrom'].drop_duplicates():
        insert_site = mapping_site_.query('chrom == @chrom').copy()
        refer_info_ = all_gene_Info_TSS.query('chrom == @chrom')
        insert_sites = np.array_split(insert_site, np.ceil(len(insert_site) / 5000))
        work_space.extend([pool.apply_async(cal_distance_groupby_chromosome, args=(i, refer_info_, 50000)) for i in insert_sites])
    pool.close()
    pool.join()
    mapping_site_TSS_distance = pd.concat([i.get() for i in work_space])
    mapping_site_TSS_distance['TSS_distance_class'] = mapping_site_TSS_distance.apply(grant_class_by_TSS_distance, axis=1)
    mapping_site_TSS_distance['relative_pos_TSS'] = mapping_site_TSS_distance['relative_pos_gene']  # copy gene strand context
    mapping_site_TSS_distance.loc[mapping_site_TSS_distance['TSS_distance_class'] == "0.0kb", 'relative_pos_TSS'] = 'inner' 

    mapping_site_TSS_distance['far_away_TSS'] = mapping_site_TSS_distance['distance'].apply(lambda x: 'No' if x <= 50000 else 'Yes')

    mapping_site = pd.merge(mapping_site, mapping_site_TSS_distance[['pos_label', 'relative_pos_TSS', 'TSS_distance_class', 'far_away_TSS']].drop_duplicates('pos_label'), left_on=['pos_label'], right_on=['pos_label'], how='left')

    mapping_site_TSS_distance.to_csv(os.path.join(outpath, 'temp', f'{outlabel}.relative_pos_TSS.csv'), sep='\t', index=False)
    mapping_site.to_csv(os.path.join(outpath, f'{outlabel}.2.csv'), sep='\t', index=False)

    chromatin_state_outfile = os.path.join(outpath, 'temp', f'{outlabel}.chrom_state.anno.csv')
    cmd = ['bedtools', 'intersect', f'-a {mapping_file}', f'-b {chrom_state}', '-wo', f'> {chromatin_state_outfile}']
    status, result = getstatusoutput(' '.join(cmd))
    if status:
        print(f'bedtools failed:\n{result}')
        exit()

    mapping_site_chrom_state = pd.read_csv(f'{chromatin_state_outfile}', sep='\t', header=None)
    chrom_state_column = len(mapping_site_chrom_state.columns) - 2
    overlap_column = len(mapping_site_chrom_state.columns) - 1
    mapping_site_chrom_state['pos_label'] = mapping_site_chrom_state.apply(lambda x: f'{x[0]}_{x[1]}_{x[2]}', axis=1)
    mapping_site_chrom_state_select = mapping_site_chrom_state.sort_values(['pos_label', overlap_column], ascending=[False, False]).drop_duplicates('pos_label')  # keep the chromatin state with largest overlap
    mapping_site_chrom_state_select = mapping_site_chrom_state_select.rename({chrom_state_column: 'chromatin_state'}, axis=1)
    mapping_site = pd.merge(mapping_site, mapping_site_chrom_state_select[['pos_label', 'chromatin_state']].drop_duplicates('pos_label'), left_on=['pos_label'], right_on=['pos_label'], how='left')
    mapping_site.to_csv(os.path.join(outpath, f'{outlabel}.3.csv'), sep='\t', index=False)

    cancer_gene = pd.read_csv(cancer_genelist, sep='\t')
    cancer_gene.columns = [i.replace(' ','_') for i in cancer_gene.columns]
    cancer_gene.loc[cancer_gene['Hugo_Symbol']=="NHERF1",'Hugo_Symbol'] = 'SLC9A3R1'    # alias override
    cancer_gene.loc[cancer_gene['Hugo_Symbol']=="U2AF1",'Hugo_Symbol'] = 'LOC102724594' # alias override

    cancer_genes = set(list(cancer_gene['Hugo_Symbol']))
    oncogene = cancer_gene.query('Is_Oncogene=="Yes"')['Hugo_Symbol'].to_list()
    suppressor = cancer_gene.query('Is_Tumor_Suppressor_Gene=="Yes"')['Hugo_Symbol'].to_list()

    work_space = []
    pool = Pool(thread)
    for chrom in mapping_site_['chrom'].drop_duplicates():
        insert_site = mapping_site_.query('chrom == @chrom').copy()
        refer_info_ = all_gene_Info_transcript.query('chrom == @chrom and gene_name in @cancer_genes')
        insert_sites = np.array_split(insert_site, np.ceil(len(insert_site) / 5000))
        work_space.extend([pool.apply_async(cal_distance_groupby_chromosome, args=(i, refer_info_, 50000)) for i in insert_sites])
    pool.close()
    pool.join()
    mapping_site_farAway_cancerGene = pd.concat([i.get() for i in work_space])
    mapping_site_farAway_cancerGene.to_csv(os.path.join(outpath, 'temp', f'{outlabel}_farAway_cancerGene.csv'), sep='\t', index=False)
    mapping_site_farAway_cancerGene['far_away_cancerGene'] = mapping_site_farAway_cancerGene['distance'].apply(lambda x: 'No' if x <= 300000 else 'Yes')  # 300kb

    mapping_site = pd.merge(mapping_site, mapping_site_farAway_cancerGene[['pos_label', 'far_away_cancerGene']].drop_duplicates('pos_label'), left_on=['pos_label'], right_on=['pos_label'], how='left')
    mapping_site.to_csv(os.path.join(outpath, f'{outlabel}.4.0.csv'), sep='\t', index=False)

    work_space = []
    pool = Pool(thread)
    for chrom in mapping_site_['chrom'].drop_duplicates():
        insert_site = mapping_site_.query('chrom == @chrom').copy()
        refer_info_ = all_gene_Info_transcript.query('chrom == @chrom and gene_name in @oncogene')
        try:
            insert_sites = np.array_split(insert_site, np.ceil(len(insert_site) / 5000))
        except Exception as error:
            print('Failed to split tasks...')
            raise RuntimeError('Failed to split insert_site chunks') from error
        work_space.extend([pool.apply_async(cal_distance_groupby_chromosome, args=(i, refer_info_, get_step(chrom,len(refer_info_)))) for i in insert_sites])
    pool.close()
    pool.join()
    mapping_site_In_oncogene = pd.concat([i.get() for i in work_space])
    mapping_site_In_oncogene.to_csv(os.path.join(outpath, 'temp', f'{outlabel}_In_oncogene.csv'), sep='\t', index=False)
    mapping_site_In_oncogene['In_oncogene'] = mapping_site_In_oncogene.apply(grant_class_by_cancerGene,axis=1)  
    mapping_site = pd.merge(mapping_site, mapping_site_In_oncogene[['pos_label', 'In_oncogene']].drop_duplicates('pos_label'), left_on=['pos_label'], right_on=['pos_label'], how='left')
    mapping_site.to_csv(os.path.join(outpath, f'{outlabel}.4.1.csv'), sep='\t', index=False)

    work_space = []
    pool = Pool(thread)
    for chrom in mapping_site_['chrom'].drop_duplicates():
        insert_site = mapping_site_.query('chrom == @chrom').copy()
        refer_info_ = all_gene_Info_transcript.query('chrom == @chrom and gene_name in @suppressor')
        insert_sites = np.array_split(insert_site, np.ceil(len(insert_site) / 5000))
        work_space.extend([pool.apply_async(cal_distance_groupby_chromosome, args=(i, refer_info_, get_step(chrom,len(refer_info_)))) for i in insert_sites])
    pool.close()
    pool.join()
    mapping_site_In_suppressor = pd.concat([i.get() for i in work_space])
    mapping_site_In_suppressor.to_csv(os.path.join(outpath, 'temp', f'{outlabel}_In_suppressor.csv'), sep='\t', index=False)
    mapping_site_In_suppressor['In_suppressor'] = mapping_site_In_suppressor.apply(grant_class_by_cancerGene,axis=1)  
    mapping_site = pd.merge(mapping_site, mapping_site_In_suppressor[['pos_label', 'In_suppressor']].drop_duplicates('pos_label'), left_on=['pos_label'], right_on=['pos_label'], how='left')
    mapping_site.to_csv(os.path.join(outpath, f'{outlabel}.4.2.csv'), sep='\t', index=False)

    miRNA = pd.read_csv(miRNA_file, sep='\t', header=None)
    miRNA.columns = ['chrom', 'start', 'end', 'gene_name', 'category']
    miRNA['strand'] = '+' 

    work_space = []
    pool = Pool(thread)
    for chrom in mapping_site_['chrom'].drop_duplicates():
        insert_site = mapping_site_.query('chrom == @chrom')
        refer_info_ = miRNA.query('chrom == @chrom and category=="miRNA"').copy()
        insert_sites = np.array_split(insert_site, np.ceil(len(insert_site) / 5000))
        work_space.extend([pool.apply_async(cal_distance_groupby_chromosome, args=(i, refer_info_, 50000)) for i in insert_sites])
    pool.close()
    pool.join()

    mapping_site_farAway_miRNA = pd.concat([i.get() for i in work_space])
    mapping_site_farAway_miRNA['far_away_miRNA'] = mapping_site_farAway_miRNA['distance'].apply(lambda x: 'No' if x <= 300000 else 'Yes')  # 300kb

    mapping_site = pd.merge(mapping_site, mapping_site_farAway_miRNA[['pos_label', 'far_away_miRNA']].drop_duplicates('pos_label'), left_on=['pos_label'], right_on=['pos_label'], how='left')
    mapping_site.to_csv(os.path.join(outpath, f'{outlabel}.5.csv'), sep='\t', index=False)

    conserved_region_outfile = os.path.join(outpath, 'temp', f'{outlabel}.conserved_region.anno.csv')
    cmd = ['bedtools', 'intersect', f'-a {mapping_file}', f'-b {conserved_region_file}', '-wo', f'> {conserved_region_outfile}']
    status, result = getstatusoutput(' '.join(cmd))
    if status:
        print(f'bedtools failed:\n{result}')
        exit()

    mapping_site_conserved_region = pd.read_csv(f'{conserved_region_outfile}', sep='\t', header=None)
    mapping_site_conserved_region['pos_label'] = mapping_site_conserved_region.apply(lambda x: f'{x[0]}_{x[1]}_{x[2]}', axis=1)
    mapping_site.loc[~mapping_site['pos_label'].isin(mapping_site_conserved_region['pos_label']), 'out_conserved_region'] = 'Yes'
    mapping_site.loc[mapping_site['pos_label'].isin(mapping_site_conserved_region['pos_label']), 'out_conserved_region'] = 'No'
    mapping_site.to_csv(os.path.join(outpath, f'{outlabel}.6.csv'), sep='\t', index=False)

    print(f'running time:{(time.time() - start_time)/60 :.2f} minutes.')
