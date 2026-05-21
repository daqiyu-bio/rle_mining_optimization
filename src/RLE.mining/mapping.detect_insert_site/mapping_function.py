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
#   Helper functions for insertion-site detection and analysis.                             #
# ------------------------------------------------------------------------------------------#

import re
import os
import pysam
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Seq, SeqIO
from function import read_psl, run_shell, read_fasta
from blat_run_analysis import blat_match_block_indels


def get_snapgene_seq(input_file):
    try:
        with open(input_file, 'rb') as dna_file:
            sequences = SeqIO.parse(dna_file, 'snapgene')
            seq = [i.seq.upper() for i in sequences][0]
    except:
        with open(input_file, 'r') as dna_file:
            sequences = SeqIO.parse(dna_file, 'gb')
            seq = [i.seq.upper() for i in sequences][0]
    return seq


def compare_plasmid_to_genome(blat, genome_fa, plasmid_fa, outpath, outlabel):
    """
    When using a plasmid as the donor, compare the plasmid sequence with 
    the genome to identify regions of homology.
    Homologous sequences between the plasmid and genome can complicate or 
    confound the detection of true insertion events in genome.
    """
    outfile = f'{outpath}/{outlabel}.plasmid_mapping2genome.csv'
    if os.path.isfile(outfile):
        return outfile

    cmd = f'{blat} {genome_fa} {plasmid_fa} -tileSize=12 -minIdentity=95 -maxIntron=10 {outpath}/{outlabel}.plasmid_mapping2genome.psl' 
    if os.path.isfile(f'{outpath}/{outlabel}.plasmid_mapping2genome.psl'):
        pass 
    else:
        run_shell(cmd)
    match, gaps, dels = blat_match_block_indels(plasmid_fa, genome_fa, outpath, run_id=f'{outlabel}.plasmid_mapping2genome')
    match = match.reset_index(drop=True)
    reference_fa = read_fasta(genome_fa)
    for i in match.index:
        chrom = match.loc[i, 'target']
        start = match.loc[i, 'target_start']
        end = match.loc[i, 'target_end']
        strand = match.loc[i, 'strand']
        if strand == "+":
            seq = str(reference_fa[chrom].seq)[start:end]
        else:
            seq = reference_fa[chrom].seq[start:end]
            seq = str(seq.reverse_complement())
        match.loc[i, 'target_seq'] = seq.upper()


    need = ['ID', 'query', 'query_start', 'query_end', 'target', 'target_start', 'target_end', 'strand',
            'match_ratio', 'identity', 'is_full', 'is_half', 'query_size']
    if match.empty:
        match_sele = pd.DataFrame(columns=need)
    else:
        match_sele = match[need]
        plasmid_len = match['query_size'].iloc[0]
        pos_stat = {}
        for i in match.index:
            start = match.loc[i, 'query_start']
            end = match.loc[i, 'query_end']
            for j in range(start, end):
                if j in pos_stat:
                    pos_stat[j] += 1
                else:
                    pos_stat[j] = 1
        fig, axes = plt.subplots(1, 1, figsize=(16, 4))
        for i in pos_stat.keys():
            axes.vlines(i, 0, pos_stat[i])
        axes.set_xlim(0, plasmid_len)
        axes.set_ylim(0)
        unit = int(10 ** np.floor(np.log10(max(pos_stat.values()))))
        axes.set_yticks(np.array(range(max(pos_stat.values()) // unit + 2)) * unit)
        axes.set_ylabel('mapping count', size=13)
        axes.set_xlabel('position of plasmid', size=13)
        axes.set_title(f'Mapping Plasmid to Genome({outlabel})', size=15)
        plt.savefig(f'{outpath}/{outlabel}.plasmid_mapping2genome.png', dpi=300, facecolor='white', bbox_inches='tight')
        plt.close(fig)
    match_sele.to_csv(outfile, sep='\t', index=False)
    return outfile


def build_bwa_index(bwa, genome_fa, plasmid_fa, outpath, outlabel):
    merge_fa = os.path.join(outpath, outlabel + '.ref.fa')
    check_file = os.path.join(outpath, 'check_file.txt')
    if os.path.isfile(check_file):
        return merge_fa
    else:
        cmd = f'cat {genome_fa} {plasmid_fa} > {merge_fa}'
        run_shell(cmd)
        cmd = f'{bwa} index -p {merge_fa} {merge_fa} && md5sum {merge_fa}.* > {check_file}'
        run_shell(cmd)
        return merge_fa



def NGS_data_qc(umi_tools,cutadapt, fastp, R1, R2, outpath, outlabel,UMI=0, ME_seq='AGATGTGTATAAGAGACAG',sequenator='Illumina'):
    """
    1. Remove the ME sequence introduced during Tn5 library preparation. ME typically
        appears at the 5' end of read1, but when the fragment is short and sequencing
        reads through, a reverse complementary sequence of ME may also appear at the 3' end of read2.
    2. Perform quality control, calculating the total data amount and the proportion of Q30 bases.
    """
    if UMI > 0:
        # detect UMI and ME sequence
        cmd = [f'{umi_tools} extract --extract-method=regex --bc-pattern="^(?P<umi_seq>.{{{UMI}}})(?P<ME_seq>{ME_seq})"']
        cmd += ['--quality-filter-threshold=20 --quality-encoding=phred33']
        if sequenator in ['MGI']:
            cmd += ['--ignore-read-pair-suffixes'] # BGI (MGI) and Illumina sequencing instruments have differences in read ID naming conventions.
       
        else:
            pass
        cmd += [f'-I {R1} -S {outpath}/{outlabel}.UMI.R1.fastq.gz']
        cmd += [f'--read2-in {R2} --read2-out {outpath}/{outlabel}.UMI.R2.fastq.gz']
        cmd += [f'--filtered-out {outpath}/{outlabel}.no_UMI.R1.fastq.gz --filtered-out2 {outpath}/{outlabel}.no_UMI.R2.fastq.gz']
        cmd += [f' > {outpath}/{outlabel}.UMI.log 2>&1']
        run_shell(' '.join(cmd))
        R1 = f'{outpath}/{outlabel}.UMI.R1.fastq.gz'
        R2 = f'{outpath}/{outlabel}.UMI.R2.fastq.gz'
        print(f'The UMI sequence({UMI}nt) in {R1} {R2} have been extracted and moved to readID.')

    ME_seq_reverse = str(Seq.Seq(ME_seq).reverse_complement())
    cmd1 = f'{cutadapt} -e 2 -g {ME_seq} -o {outpath}/{outlabel}.delME.R1.fastq.gz {R1} > {outpath}/{outlabel}.delME.R1.fastq.gz.log'
    cmd2 = f'{cutadapt} -e 2 -a {ME_seq_reverse} -o {outpath}/{outlabel}.delME.R2.fastq.gz {R2} > {outpath}/{outlabel}.delME.R2.fastq.gz.log'
    run_shell(cmd1)
    run_shell(cmd2)

    cmd3 = [f'{fastp} -i {outpath}/{outlabel}.delME.R1.fastq.gz -I {outpath}/{outlabel}.delME.R2.fastq.gz']
    cmd3 += [f'-o {outpath}/{outlabel}.delME.R1.qc.fastq.gz -O {outpath}/{outlabel}.delME.R2.qc.fastq.gz']
    cmd3 += ['-g']
    cmd3 += [f'-h {outpath}/{outlabel}.qc.html -j {outpath}/{outlabel}.qc.json']
    run_shell(' '.join(cmd3))
    return f'{outpath}/{outlabel}.delME.R1.qc.fastq.gz', f'{outpath}/{outlabel}.delME.R2.qc.fastq.gz', f'{outpath}/{outlabel}.qc.json'


def transfer_UMI2tag(bam,outpath,outlabel):
    """
    Extract the UMI from the read ID and add the UMI to the SAM file.
    """
    outbam = os.path.join(outpath, outlabel + '.RXtag.bam')
    bamR = pysam.AlignmentFile(bam, 'rb')
    bamW = pysam.AlignmentFile(outbam, 'wb', template=bamR)

    for read in bamR:
        tmp = read.qname.split('_')
        umi = tmp[-1]
        read.set_tag('RX', umi)
        read.qname = '_'.join(tmp[:-1])
        bamW.write(read)
    bamR.close()
    bamW.close()
    return outbam


def bwa_mem_align(bwa, java, picard, samtools, genome_fa, R1, R2, outpath, outlabel,UMI=0,min_score=20):
    """
    Aligns sequencing reads to both the reference genome and the donor sequence, then processes the resulting BAM files through steps such as sorting, duplicate marking, data cleaning, and statistics generation. 
    Supports sequencing data that includes Unique Molecular Identifiers (UMIs).
    By default, sets the minimum alignment score to 20, which is compatible with primers located very close to the ends of transposition templates.
    """

    # alignment
    cmd = [f'{bwa} mem -t 16 -T {min_score} -M -Y -q -R "@RG\\tID:{outlabel}\\tPL:***\\tLB:{outlabel}\\tSM:{outlabel}\\tPG:bwa"']
    cmd += [f'{genome_fa} {R1} {R2} > {outpath}/{outlabel}.sam']
    run_shell(' '.join(cmd))

    # sorting
    cmd = [f'{java} -jar {picard} SortSam --SORT_ORDER coordinate --CREATE_INDEX true  --VALIDATION_STRINGENCY SILENT']
    cmd += [f'-I {outpath}/{outlabel}.sam -O {outpath}/{outlabel}.02.sort.bam']
    run_shell(' '.join(cmd))
    
    # marking duplicates
    if UMI > 0:
        bam_with_RX = transfer_UMI2tag(f'{outpath}/{outlabel}.02.sort.bam',outpath,outlabel)
        cmd = [f'{java} -jar {picard} MarkDuplicates --TAG_DUPLICATE_SET_MEMBERS true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT']
        cmd += [f'-I {bam_with_RX} -O {outpath}/{outlabel}.03.markdup_picard.bam -M {outpath}/{outlabel}.02.markdup_picard.metrics']
        cmd += ['--MOLECULAR_IDENTIFIER_TAG RX']
    else:
        cmd = [f'{java} -jar {picard} MarkDuplicates --TAG_DUPLICATE_SET_MEMBERS true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT']
        cmd += [f'-I {outpath}/{outlabel}.02.sort.bam -O {outpath}/{outlabel}.03.markdup_picard.bam -M {outpath}/{outlabel}.02.markdup_picard.metrics']
    run_shell(' '.join(cmd))

    # cleaning:
    # 4 read unmapped ; 
    # 8 mate unmapped ; 
    # 256  not primary alignment ; 
    # 512  read fails platform/vendor quality checks ; 
    # 1024 read is PCR or optical duplicate
    cmd = f'{samtools} view -bS -F 4 -F 8 -F 256 -F 512 -F 1024 {outpath}/{outlabel}.03.markdup_picard.bam > {outpath}/{outlabel}.04.sort.rmdup.bam'
    run_shell(cmd)

    run_shell(f'{samtools} index {outpath}/{outlabel}.04.sort.rmdup.bam')
    return f'{outpath}/{outlabel}.04.sort.rmdup.bam', f'{outpath}/{outlabel}.03.markdup_picard.bam'


def get_DS_tag(read):
    try:
        DS = int(read.get_tag('DS'))
    except:
        DS = 1 
    return DS


def get_DI_tag(read):
    try:
        DI = int(read.get_tag('DI'))
    except:
        DI = read.qname
    return DI


def generate_aligned_pairs(pos, cigar):
    """
    Extract aligned pairs (read position vs. reference position) based on
    alignment coordinates and the CIGAR string.

    For primary alignments, these pairs are stored in the pysam AlignedSegment object. 
    However, for supplementary alignments (SA tag), the record has been removed from the original BAM file,
    so the aligned pairs should be computed based on the SA cigar and coordinates.
    """
    cigars = re.findall(r'\d+[A-Z]', cigar)
    start = 0 
    aligned_pairs = []
    for c in cigars:
        if c[-1] == 'M':
            num = int(c[:-1])
            aligned_pairs.extend([(start + i, pos + i) for i in range(num)])
            start += num
            pos += num
        elif c[-1] == 'D':
            num = int(c[:-1])
            aligned_pairs.extend([(None, pos + i) for i in range(num)])
            pos += num
        elif c[-1] in ['I', 'S', 'H']:
            num = int(c[:-1])
            aligned_pairs.extend([(start + i, None) for i in range(num)])
            start += num
        else:
            raise ValueError(f'Special cigar:{cigar},{c}')
    return aligned_pairs


def get_seq_by_cigar(seq, pos, cigar, return_orient='forward'):
    aligned_pairs = generate_aligned_pairs(pos, cigar)
    aligned_pos = [i[0] for i in aligned_pairs if i[1] != None]  

    if len(aligned_pos) > 1:  
        start = aligned_pos[0]
        end = aligned_pos[-1]
        aligned_seq = seq[start:end + 1]
    else:
        start = end = aligned_pos[0]
        aligned_seq = seq[start] 

    if return_orient == 'forward':
        return [start, end, aligned_seq, seq, aligned_pairs]
    else:
        return reverse_orient(start, end, seq, aligned_pairs)


def reverse_orient(start, end, seq, aligned_pairs):
    reverse_ = str(Seq.Seq.reverse_complement(Seq.Seq(seq)))
    start_ = len(seq) - end - 1  # 0-base
    end_ = len(seq) - start - 1  # 0-base
    new_aligned_pairs = []
    for i, j in zip(*zip(*aligned_pairs)):
        if i == None:
            new_aligned_pairs.append((None, j))
        else:
            new_aligned_pairs.append((len(seq) - i - 1, j))
    return [start_, end_, reverse_[start_:end_ + 1], reverse_, new_aligned_pairs]


def where_read_from(pos_plasmid, pos_genome):
    if np.average(pos_plasmid) <= np.average(pos_genome):
        return 'Right_junction'
    else:
        return 'Left_junction'


def complementary(base):
    base_map = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    return base_map[base]


def cal_overlap(a, b):
    a_start, a_end = sorted(a)
    b_start, b_end = sorted(b)

    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)

    overlap_length = max(overlap_end - overlap_start, 0)

    if overlap_length > 0:
        return overlap_length  
    else:
        if a_end <= b_start:
            return -(b_start - a_end)  
        else:
            return -(a_start - b_end)  