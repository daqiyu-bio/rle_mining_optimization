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
#   Refilter, analyze, and visualize insertion events according to the plasmid or donor     #
#   configuration used in the assay.                                                        #
# ------------------------------------------------------------------------------------------#

import os
import re
import numpy as np
import pandas as pd
from Bio import Seq
import matplotlib.pyplot as plt
from mapping_instance import event_stat

FLANK_BP = 100
JUNCTION_IDX = FLANK_BP
READ_PLASMID_FLANK = 15


def overlap_length(A, B):
    return max(0, min(A[1], B[1]) - max(A[0], B[0]))

def _aat_genomic_range(absolute_start, strand, s, e):
    if strand == 'minus':
        return f'{absolute_start - e}-{absolute_start - s}'
    return f'{absolute_start + s}-{absolute_start + e}'


def get_AT_rich(string, start, absolute_start, strand, offset):
    """
    Detect AAT-repeat microsatellites near a junction inside a short window.
    """
    range_ = ''
    seq1 = ''
    repeat = 0
    hits = []

    string = str(string).upper()
    for m in re.finditer(r'((?:AT|T)?)((?:AAT)+)', string, re.I):
        s, e = m.span()
        As = 0
        for char in string[e:e+2]:
            if char != 'A':
                break
            else:
                As += 1
        overlap_len = overlap_length([s, e+As], [start, start + offset])
        if overlap_len == 0:
            continue
        s = s + len(m.group(1))  
        seq1 = m.group(2)
        range_ = _aat_genomic_range(absolute_start, strand, s, e)
        repeat = int(np.floor((e - s) / 3))
        hits.append((range_, seq1, repeat, s, e, overlap_len))

    if not hits:
        return '', '', '', '', repeat

    range_, seq1, repeat, s, e, overlap_len = hits[np.argmax([h[2] for h in hits])] 
    if s >= start:
        seq2 = ''
        seq3 = string[s:e]
    else:
        seq2 = string[s:start]
        seq3 = string[start:e]

    return range_, seq1, seq2, seq3, repeat


def _reference_window_and_absolute_start(FASTA, chrom, genome_start, map_seq, strand, offset):
    if strand == 'minus':
        center = genome_start + len(map_seq) - offset
        temp = str(
            Seq.Seq.reverse_complement(FASTA[chrom][center - FLANK_BP : center + FLANK_BP].seq)
        ).upper()
        absolute_start = center + FLANK_BP
        return temp, absolute_start

    center = genome_start + offset
    temp = str(FASTA[chrom][center - FLANK_BP : center + FLANK_BP].seq.upper())
    absolute_start = center - FLANK_BP
    return temp, absolute_start


def _mate_strand_params(strand, absolute_start, window_len):
    mate_strand = 'minus' if strand == 'plus' else 'plus'
    junc_idx_rc = window_len - 1 - JUNCTION_IDX
    if strand == 'plus':
        abs_start_rc = absolute_start + window_len
    else:
        abs_start_rc = absolute_start - window_len
    return mate_strand, junc_idx_rc, abs_start_rc


def _choose_primary_at_rich(AT_rich, AT_rich_rev, strand, mate_strand):
    if len(AT_rich[1]) >= len(AT_rich_rev[1]):
        return AT_rich, AT_rich_rev, strand
    return AT_rich_rev, AT_rich, mate_strand

def get_AT_rich_seq(FASTA, SA_TE, uniq_fragment_total, fragment_total, outpath, outlabel):
    """
    Get and summarize AAT-repeat features for SA-supported genome insertions.
    """
    for index in SA_TE.index:
        chrom = SA_TE.loc[index, 'genome_name_y']
        genome_start = int(SA_TE.loc[index, 'genome_start_y'])
        strand = SA_TE.loc[index, 'insert_genome_strand_y']
        map_seq = SA_TE.loc[index, 'read_genome_y']
        overlap_len = int(SA_TE.loc[index, 'align_overlap_label_y'])
        offset, offset2 = 0, max(overlap_len,1)
        temp, absolute_start = _reference_window_and_absolute_start(
            FASTA, chrom, genome_start, map_seq, strand, offset
        )

        AT_rich = get_AT_rich(temp, JUNCTION_IDX, absolute_start, strand, offset2)

        mate_strand, junc_idx_rc, abs_start_rc = _mate_strand_params(
            strand, absolute_start, len(temp)
        )
        AT_rich_rev = get_AT_rich(
            str(Seq.Seq(temp).reverse_complement()),
            junc_idx_rc,
            abs_start_rc,
            mate_strand,
            offset2,
        )

        AT_rich_final, AT_rich_supp_info, AT_rich_strand = _choose_primary_at_rich(
            AT_rich, AT_rich_rev, strand, mate_strand
        )

        SA_TE.loc[index, 'AT_rich_seq'] = AT_rich_final[1]
        SA_TE.loc[index, 'AAT_count'] = AT_rich_final[4]
        SA_TE.loc[index, 'AT_rich_seq_L'] = AT_rich_final[2]
        SA_TE.loc[index, 'AT_rich_seq_R'] = AT_rich_final[3]
        SA_TE.loc[index, 'AT_rich_pos'] = f'{chrom}:{AT_rich_final[0]}'
        SA_TE.loc[index, 'AT_rich_strand'] = AT_rich_strand
        SA_TE.loc[index, 'AT_rich_info_supp'] = (
            f'{chrom}:{AT_rich_supp_info[0]}:{AT_rich_supp_info[1]}'
            if AT_rich_supp_info[1] != ''
            else ''
        )

        if AT_rich_final[0] == '':
            transposition_site = SA_TE.loc[index, 'insert_label_y'].split('_')
            SA_TE.loc[index, 'transposition_site'] = (
                f'{transposition_site[0]}:{transposition_site[1]}-{transposition_site[2]}'
            )
        else:
            SA_TE.loc[index, 'transposition_site'] = SA_TE.loc[index, 'AT_rich_pos']

        offset2_read = READ_PLASMID_FLANK
        start_read = eval(SA_TE.loc[index, 'read_plasmid_map_y'])[1] - offset2_read
        AT_rich_read = get_AT_rich(
            SA_TE.loc[index, 'read_seq_adj_y'],
            start_read,
            start_read,
            'plus',
            offset2_read,
        )
        SA_TE.loc[index, 'AT_rich_seq_in_read'] = AT_rich_read[1]
        SA_TE.loc[index, 'AT_rich_seq_in_read_countAAT'] = AT_rich_read[4]

    SA_TE.to_csv(os.path.join(outpath, f'{outlabel}.AT_rich_info.csv'), sep='\t', index=False)

    detail_df = SA_TE[['transposition_site', 'AT_rich_seq', 'DS_tag_y', 'AAT_count']].copy()
    detail_df[['chrom', 'pos']] = detail_df['transposition_site'].str.split(':', expand=True)
    detail_df[['start', 'end']] = detail_df['pos'].str.split('-', expand=True)
    detail_df = detail_df.drop(columns=['pos'])
    detail_df = detail_df.rename(columns={'DS_tag_y': 'all_fragment_count'})
    detail_df['uniq_fragment_count'] = 1

    detail_df = detail_df.groupby(
        ['chrom', 'start', 'end', 'AT_rich_seq', 'AAT_count']
    ).agg({
        'uniq_fragment_count': 'sum',
        'all_fragment_count': 'sum',
    }).reset_index()
    detail_df['uniq_fragment_ratio'] = detail_df['uniq_fragment_count'] / uniq_fragment_total
    detail_df['all_fragment_ratio'] = detail_df['all_fragment_count'] / fragment_total

    detail_df.to_csv(
        os.path.join(outpath, f'{outlabel}.AAT_rich_detail_uniq.csv'), index=False, sep='\t'
    )
    detail_df[['chrom', 'start', 'end', 'AT_rich_seq']].to_csv(
        os.path.join(outpath, f'{outlabel}.AAT_rich_detail_uniq.bed'),
        index=False,
        sep='\t',
        header=None,
    )

    aat_rich_uniq_fragment_ratio = SA_TE.groupby('AAT_count').apply(
        lambda x: round(x.shape[0] / SA_TE.shape[0] * 100, 2)
    ).sort_index()
    aat_rich_all_fragment_ratio = SA_TE.groupby('AAT_count').apply(
        lambda x: round(x['DS_tag_y'].sum() / SA_TE['DS_tag_y'].sum() * 100, 2)
    ).sort_index()

    ratio_df = pd.DataFrame({
        'AAT_count': aat_rich_uniq_fragment_ratio.index,
        'Uniq_Fragment_Ratio': aat_rich_uniq_fragment_ratio.values,
        'All_Fragment_Ratio': aat_rich_all_fragment_ratio.values,
    })
    ratio_df.to_csv(os.path.join(outpath, f'{outlabel}.AAT_rich_ratio.csv'), index=False, sep='\t')

    bar_width = 0.4
    fig, axes = plt.subplots(1, 1, figsize=(12, 6))
    axes.bar(
        aat_rich_all_fragment_ratio.index - bar_width / 1.9,
        aat_rich_all_fragment_ratio.to_list(),
        width=bar_width,
        label='Ratio by All Fragment',
    )
    axes.bar(
        aat_rich_uniq_fragment_ratio.index + bar_width / 1.9,
        aat_rich_uniq_fragment_ratio.to_list(),
        width=bar_width,
        label='Ratio by Uniq Fragment',
    )
    axes.set_ylabel('Ratio(%)', size=14)
    axes.set_title(
        f'Profile of {outlabel} insertion preference(total: {SA_TE["DS_tag_y"].sum():,.0f} Fragments)',
        size=14,
    )
    axes.set_xticks(aat_rich_all_fragment_ratio.index, size=14)
    axes.legend()
    axes.set_xlabel('AAT repeat number', size=14)
    plt.savefig(
        os.path.join(outpath, f'{outlabel}.AAT_rich_profile.png'),
        facecolor='white',
        bbox_inches='tight',
        dpi=300,
    )
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    from function import read_fasta

    parser = argparse.ArgumentParser(description='Analyze AT-rich sequences in genome editing data.')
    parser.add_argument('--outpath', type=str, required=True, help='Output directory path')
    parser.add_argument('--outlabel', type=str, required=True, help='Output file prefix label')
    parser.add_argument(
        '--fasta',
        type=str,
        required=True,
        help='Path to reference genome FASTA file',
        dest='fasta_file',
    )
    parser.add_argument(
        '--fragment_file',
        type=str,
        required=True,
        help='Path to fragment parse CSV file(*.Fragment_parse.csv)',
        dest='fragment_file',
    )

    args = parser.parse_args()
    outpath = args.outpath
    outlabel = args.outlabel
    fragment_file = args.fragment_file
    FASTA = read_fasta(args.fasta_file)
    TE = pd.read_csv(fragment_file, sep='\t')


    # According our assay design, we used the specific parameters to refilter R4 mapping result.
    # For other assay, you should modify the parameters according to your assay design.
    
    cutoff = -17
    index = TE.query('status_note.str.contains("trans_genome") and exceed_boundary_y<@cutoff').index
    TE.loc[index, 'status_note'] = 'genome_plasmid_recom(1)'
    index = TE.query('status_note.str.contains("trans_genome") and genome_start_y!=genome_start_y').index
    TE.loc[index, 'status_note'] = 'genome_plasmid_recom(1)'
    index = TE.query('status_note == "trans_genome(PA)"').index
    TE.loc[index, 'status_note'] = 'genome_plasmid_recom(1)'

    TE_file = f'{outpath}/{outlabel}.Fragment_parse.csv'
    TE.to_csv(TE_file, sep='\t', index=False)
    event_stat(TE_file, outpath, outlabel)

    SA_TE = TE.query('status_note.str.contains("trans_genome\(SA")').copy()
    uniq_fragment_total = TE.shape[0]
    fragment_total = TE['DS_tag_x'].sum()
    get_AT_rich_seq(FASTA, SA_TE, uniq_fragment_total, fragment_total, outpath, outlabel)
