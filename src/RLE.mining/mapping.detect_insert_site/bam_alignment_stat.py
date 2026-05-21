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
#   Summarize BAM alignment categories (mapped / duplicate / unmapped) and plot a pie chart.#
# ------------------------------------------------------------------------------------------#

import os
import json
import pysam
from sys import argv
import matplotlib.pyplot as plt
from argparse import ArgumentParser


def bam_alignment_stat(bam,outpath,outlabel):
    outnote = os.path.join(outpath, outlabel + '.bam_stat.json')
    outfig = os.path.join(outpath, outlabel + '.bam_stat.png')
    bamR = pysam.AlignmentFile(bam, 'rb')
    map_reads = []
    uniq_reads = []
    bad_reads = []
    mate_unmapped = []
    mate_unmapped_uniq = []
    for read in bamR:
        if read.is_secondary or read.is_supplementary:  
            continue

        if read.is_unmapped or read.is_qcfail:         
            bad_reads.append(read.qname)
            continue

        if read.mate_is_unmapped:                       
            mate_unmapped.append(read.qname)
            if not read.is_duplicate:
                mate_unmapped_uniq.append(read.qname) 
            continue

        map_reads.append(read.qname)                    
        if not read.is_duplicate:
            uniq_reads.append(read.qname)

    fail_reads_num = len(bad_reads)
    fail_frag_num = len(set(bad_reads) - set(map_reads))

    uniq_reads_num = len(uniq_reads)
    uniq_frag_num = len(set(uniq_reads))

    map_reads_num = len(map_reads)
    map_frag_num = len(set(map_reads))
    dup_frag_num = map_frag_num - uniq_frag_num

    Info = {'Note1': 'Here, "mapped" refers to fragments where both read1 and read2 are aligned to the reference.',
            'Note2': '"Mate unmapped fragment" refers to fragments where only one read (read1 or read2) is mapped.',
            'total fragment': map_frag_num + fail_frag_num,
            'mapped reads': map_reads_num,
            'mapped fragment': map_frag_num,
            'mapped reads(Fragment dedup)': uniq_reads_num,
            'mapped fragment(Fragment dedup)': uniq_frag_num,
            'dup fragment': dup_frag_num,
            'fail reads': fail_reads_num,
            'fail fragment': fail_frag_num,
            'mate unmapped fragment': len(mate_unmapped),
            'mate unmapped fragment(Fragment dedup)': len(mate_unmapped_uniq),
            'bam': bam
            } 
    with open(outnote, 'w', encoding='utf-8') as f:
        json.dump(Info, f, ensure_ascii=False, indent=4)

    unmaped = fail_frag_num + len(mate_unmapped_uniq)
    mapped_uniq = uniq_frag_num
    mapped_dup = dup_frag_num
    fig, axes = plt.subplots(1, 1, figsize=(5, 5), facecolor='white')
    axes.pie([unmaped, mapped_uniq, mapped_dup], explode=[0, 0.1, 0], colors=['red', '#27AE60', '#3498DB'],
             autopct='%.2f%%',
             labels=['unmapped', 'mapped(uniq)', 'mapped(dup)'], textprops={'fontsize': 14})
    axes.set_title(f'Alignment Stat({outlabel})', y=0.94, size=20)
    plt.savefig(outfig, dpi=300, bbox_inches='tight', facecolor='white')
    return outnote

if __name__ == "__main__":
    parser = ArgumentParser(description='Summarize reads in a BAM file by alignment category (mapped, duplicate, unmapped).')
    parser.add_argument('--bam', help='bam file.', required=True)
    parser.add_argument('--outpath', help='specify the output directory.', required=True)
    parser.add_argument('--outlabel', help='specify the output file label.', required=True)
    params = parser.parse_args()
    bam = params.bam
    outpath = params.outpath
    outlabel = params.outlabel
    outfile = bam_alignment_stat(bam, outpath, outlabel)
