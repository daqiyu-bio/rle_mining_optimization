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
#   Read-pair parsing, insertion-site detection, fragment classification, and reporting     #
#   for transposon mapping.                                                                 #
# ------------------------------------------------------------------------------------------#

import os
import re
import copy
import json
import tqdm
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import Seq
from function import read_fasta
from mapping_function import get_DS_tag, get_DI_tag, get_seq_by_cigar, where_read_from
from mapping_function import complementary, cal_overlap
import line_profiler


def flank_num(number,flank=10):
    number = int(number)
    return [max(0,number-flank),number + flank]

@line_profiler.profile
def detect_insert_site(read, donor_coordinate, FASTA, TSD='', TSD_len=0):
    """
    Parse aligned read for transposon mapping: infer plasmid/genome origin and candidate insertion features.
    """

    if len(TSD) > 0:
        TSD_len = len(TSD)
        TSD_mode = 'fixed_seq'  
    else:
        TSD_mode = 'random'     

    donor_coordinate = copy.copy(donor_coordinate)
    for donor_type in donor_coordinate: 
        for pos in donor_coordinate[donor_type]:
            donor_coordinate[donor_type][pos] = [donor_coordinate[donor_type][pos][0], donor_coordinate[donor_type][pos][-1]]
    
    ALIGN_RATIO_THRESHOLD = 0.60

    chrom = read.reference_name
    start = read.reference_start  
    end = read.reference_end      
    strand = '-' if read.is_reverse else '+'
    readN = 'read1' if read.is_read1 else 'read2'
    cigar = read.cigarstring
    mapq = read.mapq
    readQ = np.average(read.query_qualities)
    read_len = len(read.query_qualities)
    if read.is_reverse:
        seq = str(Seq.Seq.reverse_complement(Seq.Seq(read.seq)))
        seq_reverse = read.seq
        seq_in_bam = seq_reverse
    else:
        seq = read.seq
        seq_reverse = str(Seq.Seq.reverse_complement(Seq.Seq(read.seq)))
        seq_in_bam = seq
    align_ratio = (end - start) / len(seq)    
    content = {         
        'read_ID': read.qname,
        'read_flag': read.flag,
        'readN': readN,
        'read_seq_adj': seq_in_bam,
        'read_len': read_len,
        'readQ': readQ,
        'mapq': mapq,
        'cigar': cigar,
        'chrom': chrom,
        'map_start': start,
        'map_end': end,
        'strand': strand,
        'DI_tag': get_DI_tag(read),
        'DS_tag': get_DS_tag(read), }

    if read.has_tag('SA'):
        SAs = read.get_tag('SA').split(';')
    else:
        if chrom in donor_coordinate:   
            content['read_source'] = 'Single(plasmid)'  
            if np.ptp(donor_coordinate[chrom]['Right']) > 0:
                d1 = cal_overlap(flank_num(donor_coordinate[chrom]['Right'][0]), [start, end])
            else:
                d1 = 0
            if np.ptp(donor_coordinate[chrom]['Left']) > 0:
                d2 = cal_overlap(flank_num(donor_coordinate[chrom]['Left'][1]), [start, end])
            else:
                d2 = 0
            if d1 > d2 and d1 > 0:
                content['cover_primer'] = d1
                content['exceed_boundary'] = end - donor_coordinate[chrom]['Right'][1]
                content['read_source2'] = 'Right_junction'
            elif d2 > d1 and d2 > 0:
                content['cover_primer'] = d2
                content['exceed_boundary'] = donor_coordinate[chrom]['Left'][0] - start
                content['read_source2'] = 'Left_junction'
            else:
                content['cover_primer'] = 0
                content['exceed_boundary'] = 9999     # maximum value used for downstream filtering
                content['read_source2'] = ''
        else:
            content['read_source'] = 'Single(genome)' 

        if align_ratio < ALIGN_RATIO_THRESHOLD:  
            content['read_source'] = 'Multi(other)'  

        return [content]


    # detect insertion site from SA splits
    return_content = []
    for SA in SAs:
        if SA == '':
            continue
        chrom_, pos_, strand_, cigar_, mapq_, _ = SA.split(',')
        pos_ = int(pos_) - 1                              # pysam primary is 0-based; SA positions are 1-based

        if strand_ == "+":
            temp_seq = get_seq_by_cigar(seq, pos_, cigar_)
        else:
            temp_seq = get_seq_by_cigar(seq_reverse, pos_, cigar_)

        # Genome rearrangement
        if len(set(donor_coordinate) & set([chrom, chrom_])) == 0: 
            content_ = copy.deepcopy(content)
            content_['SA'] = SA
            content_['read_source'] = 'Multi(genome)'  # genome-derived SA split
            align_start, align_end = pd.DataFrame(temp_seq[4])[1].agg(['min', 'max']).astype(int).to_list()
            content_add = {
                'genome_name': chrom_,
                'genome_start': align_start,
                'genome_mapL': len(temp_seq[2]),
                'genome_align_cagar': cigar_,
                'read_genome_align': [align_start, align_end]
            }
            content_.update(content_add)
            return_content.append(content_)
            continue

        insert = {chrom + '_primary': [chrom, start, strand, cigar, mapq],          
                  chrom_: [chrom_, int(pos_), strand_, cigar_, mapq_]}              

        # Chimeric alignment between plasmid and genome/plasmid
        if chrom == chrom_: 
            check_region1 = [start, end]
            check_region2 = pd.DataFrame(temp_seq[4])[1].agg(['min', 'max']).to_list()
            record = []
            for region_pos in donor_coordinate[chrom]:
                if np.ptp(donor_coordinate[chrom][region_pos]) == 0: 
                    continue    

                if region_pos == 'Right':
                    region = flank_num(donor_coordinate[chrom][region_pos][0])
                else:
                    region = flank_num(donor_coordinate[chrom][region_pos][1])

                d1 = cal_overlap(region, check_region1)
                d2 = cal_overlap(region, check_region2)
                record.append([region_pos, 'region1', d1])
                record.append([region_pos, 'region2', d2])
            sele_plasmid = pd.DataFrame(record).sort_values(2).iloc[-1].to_dict()
            if sele_plasmid[1] == 'region1': 
                plasmid = chrom + '_primary'
                plasmid_name = chrom
                genome = chrom_
            else:                
                plasmid = chrom_
                plasmid_name = chrom_
                genome = chrom + '_primary'
        else: 
            if chrom in donor_coordinate:
                plasmid = chrom + '_primary'
                plasmid_name = chrom
                genome = chrom_
            else:
                plasmid = chrom_
                plasmid_name = chrom_   
                genome = chrom + '_primary'

        if insert[plasmid][2] == '+': 
            reference = seq
        else:
            reference = seq_reverse
        content['read_seq_adj'] = reference

        plasmid_seq = get_seq_by_cigar(reference, insert[plasmid][1], insert[plasmid][3]) 
        orient = {True: 'forward', False: 'reverse'}
        insert_strand = {True: '+', False: '-'}
        same_strand = insert[plasmid][2] == insert[genome][2]  
        if insert[genome][2] == '+':                           
            genome_seq = get_seq_by_cigar(seq, insert[genome][1], insert[genome][3], return_orient=orient[same_strand])
        else:
            genome_seq = get_seq_by_cigar(seq_reverse, insert[genome][1], insert[genome][3], return_orient=orient[same_strand])

        read_base = {i:plasmid_seq[3][i] for i in range(len(plasmid_seq[3]))}
        read_align_df = pd.DataFrame([read_base]).T
        read_align_df['align_plasmid'] = read_align_df.index.map({i:j for i,j in plasmid_seq[4]})
        read_align_df['align_genome'] = read_align_df.index.map({i:j for i,j in genome_seq[4]})
        read_align_df = read_align_df.reset_index()
        read_align_df.columns = ['pos','base','align_plasmid','align_genome']
        read_align_df = read_align_df[['base','pos','align_plasmid','align_genome']]

        read_align_df['align_plasmid'] = read_align_df['align_plasmid'].fillna(-99).astype(int)
        read_align_df['align_genome'] = read_align_df['align_genome'].fillna(-99).astype(int)

        plasmid_mask = read_align_df['align_plasmid'] != -99
        genome_mask = read_align_df['align_genome'] != -99

        # alignment on plasmid
        if plasmid_mask.any():
            plasmid_pos = read_align_df.loc[plasmid_mask, 'align_plasmid']
            plasmid_bases = np.array([FASTA[insert[plasmid][0]][p].upper() for p in plasmid_pos])
            n_mask = plasmid_bases == 'N'
            if n_mask.any():
                full_n_mask = pd.Series(False, index=read_align_df.index)
                full_n_mask.loc[plasmid_seq] = n_mask
                read_align_df.loc[full_n_mask, 'align_plasmid'] = -99
                plasmid_mask &= ~full_n_mask
            read_align_df.loc[plasmid_mask, 'plasmid_base'] = plasmid_bases[~n_mask]
            plasmid_maps = read_align_df.loc[plasmid_mask, 'align_plasmid'].tolist()

        # alignment on genome
        if genome_mask.any():
            genome_pos = read_align_df.loc[genome_mask, 'align_genome']
            genome_bases = np.array([FASTA[insert[genome][0]][p].upper() for p in genome_pos])
            if not same_strand:
                genome_bases = np.array([complementary(b) for b in genome_bases])
            n_mask = genome_bases == 'N'
            if n_mask.any():
                full_n_mask = pd.Series(False, index=read_align_df.index)
                full_n_mask.loc[genome_mask] = n_mask
                read_align_df.loc[full_n_mask, 'align_genome'] = -99
                genome_mask &= ~full_n_mask
            read_align_df.loc[genome_mask, 'genome_base'] = genome_bases[~n_mask]
            genome_maps = read_align_df.loc[genome_mask, 'align_genome'].tolist()

        filtered = read_align_df[read_align_df['align_plasmid'] != -99]['pos']
        plasmid_seq[0], plasmid_seq[1] = filtered.min(), filtered.max() + 1
        plasmid_seq[2] = ''.join(read_align_df.loc[read_align_df['align_plasmid'] != -99, 'base'].values)

        filtered = read_align_df[read_align_df['align_genome'] != -99]['pos']
        genome_seq[0], genome_seq[1] = filtered.min(),filtered.max() + 1
        genome_seq[2] = ''.join(read_align_df.loc[read_align_df['align_genome'] != -99, 'base'].values)

        read_align_df['plasmid_base'] = read_align_df['plasmid_base'].fillna('_')
        read_align_df['genome_base'] = read_align_df['genome_base'].fillna('_')
        read_align_df = read_align_df[['base', 'pos', 'align_plasmid', 'align_genome', 'plasmid_base', 'genome_base']]

        # Infer Left_junction vs Right_junction from plasmid vs genome coordinates
        read_source2 = where_read_from(plasmid_seq[:2], genome_seq[:2]) 
        if read_source2 == 'Right_junction' and np.ptp(donor_coordinate[plasmid_name]['Right']) != 0:
            cover_primer = cal_overlap(flank_num(donor_coordinate[plasmid_name]['Right'][0]), [min(plasmid_maps), max(plasmid_maps) + 1])
            exceed_boundary = max(plasmid_maps) + 1 - donor_coordinate[plasmid_name]['Right'][1]
        elif read_source2 == 'Left_junction' and np.ptp(donor_coordinate[plasmid_name]['Left']) != 0:
            cover_primer = cal_overlap(flank_num(donor_coordinate[plasmid_name]['Left'][1]), [min(plasmid_maps), max(plasmid_maps) + 1])
            exceed_boundary = donor_coordinate[plasmid_name]['Left'][0] - min(plasmid_maps)
        else:
            cover_primer = 0
            exceed_boundary = 9999      
        overlap = read_align_df.loc[(read_align_df['align_plasmid']!=-99) & (read_align_df['align_genome']!=-99)]
        pattern_seq = ''
        pattern_seq_map_genome = ''
        if not overlap.empty:
            overlap_start, overlap_end = overlap.index[[0, -1]]
            overlap_seq = ''.join(read_align_df.loc[overlap_start:overlap_end]['base'].to_list())
            status = f'overlap:{overlap_seq}'
            status_len = len(overlap)

            if TSD_mode == 'fixed_seq':  
                if read_source2 == 'Right_junction':
                    plasmid_end = donor_coordinate[plasmid_name]['Right'][-1]
                    sele = read_align_df.loc[ (plasmid_end-TSD_len<=read_align_df['align_plasmid']) & (read_align_df['align_plasmid']<plasmid_end) ]
                    if set(sele.index) <= set(overlap.index)  and set(sele.index) != set():  # set subset check
                        pattern_seq = ''.join(sele['base'].to_list())
                        pattern_seq_map_genome = ''.join(sele['genome_base'].to_list())
                        pattern_seq_map_genome_start = sele['align_genome'].min()
                        pattern_seq_map_genome_end = sele['align_genome'].max() + 1
                        insert_label = f'{insert[genome][0]}_{pattern_seq_map_genome_start}_{pattern_seq_map_genome_end}_{insert_strand[same_strand]}'
                    else:
                        pattern_seq_map_genome_start = overlap['align_genome'].min()
                        pattern_seq_map_genome_end = overlap['align_genome'].max() + 1
                        insert_label = f'{insert[genome][0]}_{pattern_seq_map_genome_start}_{pattern_seq_map_genome_end}_{insert_strand[same_strand]}'
                else:
                    plasmid_start = donor_coordinate[plasmid_name]['Left'][0]
                    sele = read_align_df.loc[ (plasmid_start<=read_align_df['align_plasmid']) & (read_align_df['align_plasmid']<plasmid_start+TSD_len)]
                    if set(sele.index) <= set(overlap.index)  and set(sele.index) != set():
                        pattern_seq = ''.join(sele['base'].to_list())
                        pattern_seq_map_genome = ''.join(sele['genome_base'].to_list())
                        pattern_seq_map_genome_start = sele['align_genome'].min()
                        pattern_seq_map_genome_end = sele['align_genome'].max() + 1
                        insert_label = f'{insert[genome][0]}_{pattern_seq_map_genome_start}_{pattern_seq_map_genome_end}_{insert_strand[same_strand]}'
                    else:
                        pattern_seq_map_genome_start = overlap['align_genome'].min()
                        pattern_seq_map_genome_end = overlap['align_genome'].max() + 1
                        insert_label = f'{insert[genome][0]}_{pattern_seq_map_genome_start}_{pattern_seq_map_genome_end}_{insert_strand[same_strand]}'
            else:
                pattern_seq_map_genome_start = overlap['align_genome'].min()
                pattern_seq_map_genome_end = overlap['align_genome'].max() + 1
                insert_label = f'{insert[genome][0]}_{pattern_seq_map_genome_start}_{pattern_seq_map_genome_end}_{insert_strand[same_strand]}'

        else:
            if read_source2 == 'Right_junction':
                gap_start = read_align_df.query('align_plasmid!=-99').index[-1] + 1
                gap_end = read_align_df.query('align_genome!=-99').index[0] - 1
                gap_seq = ''.join(read_align_df.loc[gap_start:gap_end]['base'].to_list())
                status = f'gap:{gap_seq}'
                status_len = -len(gap_seq)
                break_point = read_align_df.loc[gap_end + 1, 'align_genome']
            else:
                gap_start = read_align_df.query('align_genome!=-99').index[-1] + 1
                gap_end = read_align_df.query('align_plasmid!=-99').index[0] - 1
                gap_seq = ''.join(read_align_df.loc[gap_start:gap_end]['base'].to_list())
                status = f'gap:{gap_seq}'
                status_len = -len(gap_seq)
                break_point = read_align_df.loc[gap_start - 1, 'align_genome']

            insert_label = f'{insert[genome][0]}_{break_point}_{break_point+1}_{insert_strand[same_strand]}'

        filtered = read_align_df[read_align_df['align_plasmid'] != -99]['align_plasmid']
        read_plasmid_align = [filtered.min(),filtered.max()+1]

        filtered = read_align_df[read_align_df['align_genome'] != -99]['align_genome']
        read_genome_align = [filtered.min(),filtered.max()+1]

        content_ = copy.deepcopy(content)
        if insert[genome][0] in donor_coordinate:
            read_source = 'Multi(plasmid)'
        else:
            read_source = 'Multi(plasmid+genome)'

        content_add = {
            'SA': SA,
            'read_source': read_source,                                   
            'read_source2': read_source2,                                 
            'plasmid_name': insert[plasmid][0],
            'plasmid_start': insert[plasmid][1],                          
            'plasmid_mapL': len(plasmid_seq[2]),
            'plasmid_align_cigar': insert[plasmid][3],
            'cover_primer': cover_primer,                                  
            'exceed_boundary': exceed_boundary,                            
            'genome_name': insert[genome][0],
            'genome_start': insert[genome][1],                            
            'genome_mapL': len(genome_seq[2]),
            'genome_align_cagar': insert[genome][3],
            'read_plasmid_map': tuple(plasmid_seq[0:2]),
            'read_plasmid': plasmid_seq[2],
            'read_genome_map': tuple(genome_seq[0:2]),
            'read_genome': genome_seq[2],
            'insert_genome_strand': 'plus' if same_strand else 'minus',
            'align_overlap': status,
            'align_overlap_label': status_len,
            'insert_concern_pattern': pattern_seq,                        
            'pattern_map2genome': pattern_seq_map_genome,
            'insert_label': insert_label,                                 
            'read_plasmid_align':read_plasmid_align,
            'read_genome_align': read_genome_align,
        }
        content_.update(content_add)
        return_content.append(content_)
    return return_content


@line_profiler.profile
def screen_insert_event(bam, donor_coordinate, reference, outpath, outlabel, TSD='', TSD_len=0):    
    """
    Scan a BAM for transposition-related events.
    """

    FASTA = read_fasta(reference)
    bamR = pysam.AlignmentFile(bam, 'rb')
    all_ret = []
    for read in tqdm.tqdm(bamR):
        ret = detect_insert_site(read, donor_coordinate, FASTA, TSD, TSD_len)
        all_ret.extend(ret)

    info = pd.DataFrame(all_ret)
    columns_order = ['read_ID', 'read_flag', 'readN', 'read_seq_adj', 'read_len', 'readQ', 'mapq', 'cigar', 'chrom',
                     'map_start', 'map_end', 'strand', 'DI_tag', 'DS_tag', 'SA', 'read_source', 'read_source2', 'plasmid_name',
                     'plasmid_start', 'plasmid_mapL', 'plasmid_align_cigar', 'cover_primer', 'exceed_boundary', 'genome_name', 'genome_start',
                     'genome_mapL', 'genome_align_cagar', 'read_plasmid_map', 'read_plasmid', 'read_genome_map', 'read_genome',
                     'insert_genome_strand', 'align_overlap', 'align_overlap_label', 'insert_concern_pattern',
                     'pattern_map2genome', 'insert_label', 'read_plasmid_align','read_genome_align']
    for col in columns_order:
        if col not in info.columns:
            info[col] = ''
    info = info[columns_order]
    info_R1 = info.query('readN=="read1"')
    info_R2 = info.query('readN=="read2"')
    Info_merge = info_R1.merge(info_R2, left_on=['read_ID'], right_on=['read_ID'], how='outer')
    outfile = os.path.join(outpath, f'{outlabel}.PE_reads_parse.csv')
    Info_merge.to_csv(outfile, sep='\t', index=False)
    return outfile


def check_hit_homo(plasmid_to_genome, chrom, start, end):
    mask = ((plasmid_to_genome['target'] == chrom) & 
            (plasmid_to_genome['target_start'] <= start) & 
            (plasmid_to_genome['target_end'] >= end))
    return not plasmid_to_genome[mask].empty


@line_profiler.profile
def assign_category(I, plasmid_to_genome):
    """
    Classify fragments into different event categories.
    """

    if I['readQ_x'] < 30 or I['readQ_y'] < 30:
        return 'low_qual'
    
    if str(I['cover_primer_y']) == 'nan' or I['cover_primer_y'] <= 0:
        return 'unexpected'
    if I['read_source_x'] in ['Multi(other)'] or I['read_source_y'] in ['Multi(other)']:
        return 'low_align_ratio'

    if I['read_source_x'] in ["Single(genome)", "Multi(genome)"] and I['read_source_y'] in ["Single(genome)", "Multi(genome)"]:
        status = 'genome' 
        return status

    if I['read_source_x'] in ["Single(plasmid)", "Multi(plasmid)"] and I['read_source_y'] in ["Multi(plasmid)"]:
        if I['exceed_boundary_y'] <= 10 and I['align_overlap_label_y'] >= 0:
            status = 'trans_plasmid(SA)'
            return status
        else:
            status = 'plasmid(1)'
            return status

    if I['read_source_x'] in ["Single(plasmid)", "Multi(plasmid)"] and I['read_source_y'] in ["Single(plasmid)"]:
        if I['exceed_boundary_y'] <= 10 and cal_overlap([I['map_start_x'], I['map_end_x']], [I['map_start_y'], I['map_end_y']]) < -800:
            status = 'trans_plasmid(PA)'
            return status
        else:
            status = 'plasmid(2)'
            return status

    if I['read_source_x'] in ["Single(genome)", "Multi(genome)"] and I['read_source_y'] in ["Multi(plasmid+genome)"]:
        concordant = I['genome_name_y'] == I['chrom_x'] and cal_overlap(eval(I['read_genome_align_y']), [int(I['map_start_x']), int(I['map_end_x'])]) > -800
        
        if I['read_source_x'] in ["Multi(genome)"]:
            concordant_ = I['genome_name_y'] == I['genome_name_x'] and cal_overlap(eval(I['read_genome_align_y']), eval(I['read_genome_align_x'])) > -800
            concordant = concordant or concordant_
        genome_chrom = I['genome_name_y']
        genome_start, genome_end = eval(I['read_genome_align_y'])
        hit_homo = check_hit_homo(plasmid_to_genome, genome_chrom, genome_start, genome_end)

        if I['exceed_boundary_y'] <= 10 and I['align_overlap_label_y'] >= 0 and not hit_homo:
            if concordant:
                status = 'trans_genome(SA)'
                return status
            else:
                status = 'genome_plasmid_recom(1)'
                return status
        else:
            status = 'genome_plasmid_recom(2)'
            return status

    if I['read_source_x'] in ["Multi(plasmid+genome)"] and I['read_source_y'] in ["Multi(plasmid+genome)"]:
        concordant = I['genome_name_x'] == I['genome_name_y'] and cal_overlap(eval(I['read_genome_align_x']), eval(I['read_genome_align_y'])) > 0
        genome_chrom = I['genome_name_y']
        genome_start, genome_end = eval(I['read_genome_align_y'])
        hit_homo = check_hit_homo(plasmid_to_genome, genome_chrom, genome_start, genome_end)
        if I['exceed_boundary_x'] <= 10 and I['align_overlap_label_x'] >= 0 and I['exceed_boundary_y'] <= 10 and I['align_overlap_label_y'] >= 0 and not hit_homo:
            if concordant:
                status = 'trans_genome(SA_SA)'
                return status
            else:
                status = 'genome_plasmid_recom(3)'
                return status
        else:
            status = 'genome_plasmid_recom(4)'
            return status

    if I['read_source_x'] in ["Single(genome)", "Multi(genome)"] and I['read_source_y'] in ["Single(plasmid)"]:
        genome_chrom = I['genome_name_x']
        genome_start = I['genome_start_x'] 
        genome_end = I['genome_start_x'] + I['genome_mapL_x'] - 1
        hit_homo = check_hit_homo(plasmid_to_genome, genome_chrom, genome_start, genome_end)

        if I['exceed_boundary_y'] <= 10 and not hit_homo:
            status = 'trans_genome(PA)'
            return status
        elif I['exceed_boundary_y'] <= 10:
            status = 'trans_plasmid(PA)'
            return status
        else:
            status = 'genome_plasmid_recom(5)'
            return status


    if I['read_source_x'] in ["Multi(plasmid+genome)"] and I['read_source_y'] in ["Single(plasmid)"]:
        genome_chrom = I['genome_name_x']
        genome_start = I['genome_start_x'] 
        genome_end = I['genome_start_x'] + I['genome_mapL_x'] - 1
        hit_homo = check_hit_homo(plasmid_to_genome, genome_chrom, genome_start, genome_end)

        if not hit_homo:
            if I['cover_primer_x'] > 0  and I['exceed_boundary_y'] <= 10:
                status = 'trans_genome(SA)'        
                return status
            else:
                status = 'genome_plasmid_recom(6)' 
                return status
        else:  
            if I['cover_primer_x'] > 0  and I['exceed_boundary_y'] <= 10:
                status = 'trans_plasmid(PA)'
                return status
            else:
                status = 'plasmid(3)'
                return status
    
    status = 'other'
    return status

def batch_assign_category(subInfo, plasmid_to_genome):
    ret = {}
    for i in subInfo.index:
        status = assign_category(subInfo.loc[i].to_dict(), plasmid_to_genome)
        ret[i] = status
    return ret

def split_task(task_index_list,thread_number):
    Load = int(np.ceil(len(task_index_list)/thread_number))
    task_list = []
    for i in range(thread_number):
        task_list.append(task_index_list[i*Load:(i+1)*Load])
    return task_list



def fragment_classification(reads_parse_file, plasmid_map2genome_file, outpath, outlabel):
    plasmid_to_genome = pd.read_csv(plasmid_map2genome_file, sep='\t')
    Info = pd.read_csv(reads_parse_file, sep='\t')

    print('Classifying fragment origins...')
    thread_number = 20
    with ThreadPoolExecutor(max_workers=thread_number) as executor:
        future_category = [executor.submit(batch_assign_category, Info.loc[task], plasmid_to_genome) for task in split_task(Info.index,thread_number*3)]
        category_dict = {}
        i = 0
        total = len(future_category)
        for future in as_completed(future_category):
            try:
                category = future.result()
                category_dict.update(category)
            except Exception as error:
                print(f"Error while classifying batch: {error}", category)
                traceback.print_exc()
                raise Exception(error)
            print('progress2:',f'{round(i/total*100,2):.2f}%',end='\r')
            i += 1
    
    Info['status_note'] = Info.index.map(category_dict)

    category_order = [          
        'trans_genome(SA_SA)',
        'trans_genome(SA)',
        'trans_genome(PA)',
        'trans_plasmid(SA)',
        'trans_plasmid(PA)',
        'plasmid(1)',
        'plasmid(2)',
        'plasmid(3)',
        'genome',
        'unexpected',
        'genome_plasmid_recom(1)',
        'genome_plasmid_recom(2)',
        'genome_plasmid_recom(3)',
        'genome_plasmid_recom(4)',
        'genome_plasmid_recom(5)',
        'genome_plasmid_recom(6)',
        'low_qual',
        'low_align_ratio',
        'other', ]
    Info['status_note'] = pd.Categorical(Info['status_note'], categories=category_order, ordered=True)
    Info_choose_event = Info.sort_values(['read_ID', 'status_note']).drop_duplicates('read_ID')
    outfile = os.path.join(outpath, f'{outlabel}.Fragment_parse.csv')
    Info_choose_event.to_csv(outfile, sep='\t', index=False)
    return outfile

def event_stat(fragment_event_file, outpath, outlabel):
    """
    Summarize and plot fragment event categories.
    """

    fragment_event = pd.read_csv(fragment_event_file,sep='\t')
    fragment_event['DS_tag'] = fragment_event['DS_tag_x'].fillna(fragment_event['DS_tag_y'])
    cate_groups = {
        'trans_genome': ['trans_genome(SA_SA)', 'trans_genome(SA)', 'trans_genome(PA)'],
        'plasmid': ['plasmid(1)', 'plasmid(2)','plasmid(3)'],
        'trans_plasmid': ['trans_plasmid(SA)', 'trans_plasmid(PA)'],
        'plasmid_genome': ['genome_plasmid_recom(1)', 'genome_plasmid_recom(2)', 'genome_plasmid_recom(3)', 
                           'genome_plasmid_recom(4)', 'genome_plasmid_recom(5)', 'genome_plasmid_recom(6)'],
        'bad_data': ['low_qual','unexpected','low_align_ratio'],
        'other': ['other']
    }

    event_stat = fragment_event.groupby('status_note').agg(
        Uniq_Fragment=('DS_tag', 'count'),  
        All_Fragment=('DS_tag', 'sum')      
    ).reset_index()
    
    total_row = pd.DataFrame({
        'status_note': ['Total'],
        'Uniq_Fragment': [event_stat['Uniq_Fragment'].sum()],
        'All_Fragment': [event_stat['All_Fragment'].sum()]
    })
    event_stat = pd.concat([event_stat, total_row], ignore_index=True)
    
    event_stat['Uniq_ratio'] = event_stat['Uniq_Fragment'] / fragment_event.shape[0] * 100
    event_stat['All_ratio'] = event_stat['All_Fragment'] / fragment_event['DS_tag'].sum() * 100
    event_stat[['Uniq_ratio', 'All_ratio']] = event_stat[['Uniq_ratio', 'All_ratio']].round(2)
    
    event_stat.sort_values('All_Fragment', ascending=False).to_csv(f'{outpath}/{outlabel}.event_stat.csv', sep='\t', index=False)
    
    content = {}
    for item in cate_groups:
        group = cate_groups[item]
        uniq = fragment_event.query('status_note.isin(@group)')['DS_tag'].shape[0]
        All = fragment_event.query('status_note.isin(@group)')['DS_tag'].sum()
        content[item] = [uniq/fragment_event.shape[0]*100,All/fragment_event['DS_tag'].sum()*100,uniq,All]
    stat_to_plot = pd.DataFrame(content).apply(lambda x:round(x,4))

    bar_width = 0.35  
    index = np.arange(len(cate_groups))  
    fig,axes = plt.subplots(1,1,figsize=(10,4))
    y_offset = max(stat_to_plot.loc[0].max(),stat_to_plot.loc[1].max()) * 0.02
    axes.bar(index - bar_width/1.9 , stat_to_plot.loc[1].to_list(), width=bar_width, label='Ratio by All Fragment')
    [axes.text(x_,y_+y_offset,round(text,2),ha='center') for x_,y_,text in zip(index - bar_width/1.9 , stat_to_plot.loc[1].to_list(),stat_to_plot.loc[1].to_list())]

    axes.bar(index + bar_width/1.9 , stat_to_plot.loc[0].to_list(), width=bar_width, label='Ratio by Uniq Fragment')
    [axes.text(x_,y_+y_offset,round(text,2),ha='center') for x_,y_,text in zip(index + bar_width/1.9 , stat_to_plot.loc[0].to_list(),stat_to_plot.loc[0].to_list())]

    axes.set_ylim(0,y_offset/0.02 *1.1)
    axes.set_ylabel('Ratio(%)',size=12)
    axes.set_title(f'Mapping Data Composition ( {outlabel} total:{int(stat_to_plot.loc[3].sum()):,} Fragments )',size=12)
    axes.set_xticks(index, cate_groups,size=10) 
    axes.legend()  
    plt.savefig(f'{outpath}/{outlabel}.Mapping_data_composition.png',facecolor='white',bbox_inches='tight',dpi=300)
    
    content2 = {}
    for item in cate_groups['bad_data']:
        group = [item]
        uniq = fragment_event.query('status_note.isin(@group)')['DS_tag'].shape[0]
        All = fragment_event.query('status_note.isin(@group)')['DS_tag'].sum()
        content2[item] = [uniq/stat_to_plot.loc[2]['bad_data']*100,All/stat_to_plot.loc[3]['bad_data']*100,uniq,All]
    stat_to_plot2 = pd.DataFrame(content2).apply(lambda x:round(x,4))
    
    bar_width = 0.3 
    index = np.arange(len(cate_groups['bad_data'])) 
    fig,axes = plt.subplots(1,1,figsize=(5,4))
    y_offset = max(stat_to_plot2.loc[0].max(),stat_to_plot2.loc[1].max()) * 0.02
    axes.bar(index - bar_width/1.9 , stat_to_plot2.loc[1].to_list(), width=bar_width, label='Ratio by All Fragment')
    [axes.text(x_,y_+y_offset,round(text,2),ha='center',size=12) for x_,y_,text in zip(index - bar_width/1.9 , stat_to_plot2.loc[1].to_list(),stat_to_plot2.loc[1].to_list())]

    axes.bar(index + bar_width/1.9 , stat_to_plot2.loc[0].to_list(), width=bar_width, label='Ratio by Uniq Fragment')
    [axes.text(x_,y_+y_offset,round(text,2),ha='center',size=12) for x_,y_,text in zip(index + bar_width/1.9 , stat_to_plot2.loc[0].to_list(),stat_to_plot2.loc[0].to_list())]

    axes.set_ylim(0,y_offset/0.02 *1.1)
    axes.set_ylabel('Ratio(%)',size=12)
    axes.set_title(f'Mapping Data(bad) Composition({outlabel})',size=12)
    axes.set_xticks(index, cate_groups['bad_data'],size=12) 
    axes.legend()  
    plt.savefig(f'{outpath}/{outlabel}.Mapping_data_composition_bad.png',facecolor='white',bbox_inches='tight',dpi=300)

def data_stat(bam_stat_file,qc_file,fragment_event_file,outpath,outlabel):
    """
    Aggregate QC/BAM stats and trigger ``event_stat`` plots/tables.
    """

    with open(bam_stat_file) as fw:
        bam_stat = json.load(fw)
    with open(qc_file) as fw:
        qc_dict = json.load(fw)
    total_base = qc_dict['summary']['before_filtering']['total_bases']
    total_Fragment = qc_dict['summary']['before_filtering']['total_reads']/2
    qc30_rate = qc_dict['summary']['before_filtering']['q30_rate']

    map_fragment = bam_stat['mapped fragment']
    map_fragment_uniq = bam_stat['mapped fragment(Fragment dedup)']
    map_ratio = round(map_fragment/total_Fragment,4)
    dup_ratio = round(1 - map_fragment_uniq/map_fragment,4)

    content = {
        'total_base': total_base,
        'total_Fragment': int(total_Fragment),
        'qc30_rate': qc30_rate,
        'map_fragment': map_fragment,
        'map_ratio': map_ratio,
        'map_fragment_uniq': map_fragment_uniq,
        'dup_ratio': dup_ratio
    }
    pd.DataFrame([content]).to_csv(f'{outpath}/{outlabel}.data_stat.csv',sep='\t',index=False)
    
    event_stat(fragment_event_file, outpath, outlabel)
    return True


def show_data_dupstatus(Fragment_status,outpath,outlabel):
    """
    Plot duplication depth for the top fragments and cumulative proportion of library mass.
    """
    show_number = 2000
    trans_color = '#03A791'
    other_color = '#AF8260'
    cumsum_color= '#102E50'
    legend_bg_color = '#E4EFE7'

    temp = pd.read_csv(Fragment_status,sep='\t')
    temp['color'] = temp['status_note'].apply(lambda x: trans_color if re.match('trans_genome',x) else other_color)
    total_fragment = temp['DS_tag_x'].sum()
    temp = temp.sort_values('DS_tag_x',ascending=False).head(show_number)
    fig,axes = plt.subplots(figsize=(16,7))
    axes.bar(range(len(temp)),np.log2(np.array(temp['DS_tag_x'].to_list())),color=temp['color'])
    axes.set_ylabel('log2(count)',size=16)

    ratio = np.array(temp['DS_tag_x'])/total_fragment
    ratio_cumsum = np.cumsum(ratio)
    axes2 = axes.twinx()
    axes2.plot(range(len(ratio_cumsum)),ratio_cumsum*100,'--',color=cumsum_color,linewidth =3,label='Cumul. Proportion')
    axes2.set_ylabel('Cumulative Proportion(%)',size=16,rotation=270,labelpad=20)
    axes2.grid()
    axes2.set_ylim(0,110)
    axes.set_xlabel(f'Uniq Fragment No.(limit={min(show_number,len(temp))})',size=16)
    axes.set_title(f'Duplication Status({outlabel})',size=18,y=1.01)
    axes.set_xlim(-show_number*0.02,show_number * 1.02)
    axes2.bar(0,0,color=trans_color,label='Trans Event')
    axes2.bar(0,0,color=other_color,label='Other Event')

    axes2.legend(loc='center right',facecolor=legend_bg_color,fontsize=14)
    plt.savefig(f'{outpath}/{outlabel}.dup.png',dpi=300,bbox_inches='tight',facecolor='white')  