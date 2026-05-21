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
#   Calculate and visualize enrichment ratios for multiple insertion-site feature           #
#   categories.                                                                             #
# ------------------------------------------------------------------------------------------#

import os
import yaml
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from matplotlib.font_manager import FontProperties
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

import glob
import matplotlib as mpl
from matplotlib import font_manager

# specify the font path and add it to the font manager
# for font in glob.glob("/path/to/fonts/CALIBRI*"):
#     font_manager.fontManager.addfont(font)

# mpl.rcParams.update({
#     'pdf.fonttype': 42,            
#     'ps.fonttype': 42,             
#     'font.family': 'sans-serif',   
#     'font.sans-serif':'Calibri'
# })

def parse_params():
    parser = ArgumentParser(description='Read categorized insertion-site tables produced by insert_site_annotation.py.')
    parser.add_argument('--mapping_result', dest='mapping_file', action='store', help='result of mapping data with category', required=True)
    parser.add_argument('--sample_label', dest='sample_label', action='store', help='specify sample label matched with mapping_result.', required=True)
    parser.add_argument('--random', dest='random', action='store', help='result of simulated data with category.', required=True)
    parser.add_argument('--outpath', dest='outpath', action='store', help='specify the output directory.', required=True)
    parser.add_argument('--outlabel', dest='outlabel', action='store', help='specify the output file label.', required=True)
    parser.add_argument('--sample_order', dest='sample_order', action='store', help='specify the sample order file', default='none')
    params = parser.parse_args()
    return params


_WEIGHT_COL_PRIORITY = ('uniq_fragment_count', 'fragment_count')


def _pick_weight_column(columns):
    for c in _WEIGHT_COL_PRIORITY:
        if c in columns:
            return c
    return None


def cal_feature_ratio(site_with_feature_mark):
    """
    Compute per-category ratios for the annotated insertion-site table.
    """
    mapping_site_feature = pd.read_csv(site_with_feature_mark, sep='\t')
    weight_col = _pick_weight_column(mapping_site_feature.columns)
    use_fragments = weight_col is not None
    if use_fragments:
        mapping_site_feature[weight_col] = pd.to_numeric(
            mapping_site_feature[weight_col], errors='coerce'
        ).fillna(0)

    def row_metric(subdf):
        if use_fragments:
            return float(subdf[weight_col].sum())
        return float(len(subdf))

    def qcount(q):
        return row_metric(mapping_site_feature.query(q))

    safe_harbor_columns = ['far_away_TSS', 'far_away_cancerGene', 'far_away_miRNA', 'out_gene', 'out_conserved_region']
    mapping_site_feature['combined_safe_harbor'] = (mapping_site_feature[safe_harbor_columns] == "Yes").apply(lambda x: np.all(x), axis=1)

    chromatin_states = [
        '1_TssA',
        '2_TssAFlnk',
        '3_TxFlnk',
        '4_Tx',
        '5_TxWk',
        '6_EnhG',
        '7_Enh',
        '8_ZNF/Rpts',
        '9_Het',
        '10_TssBiv',
        '11_BivFlnk',
        '12_EnhBiv',
        '13_ReprPC',
        '14_ReprPCWk',
        '15_Quies',
    ]

    feature_count = { 
        'up_gt_50kb': qcount('relative_pos_gene=="upstream" and gene_distance_class==">50kb"'),
        'up_lt_50kb': qcount('relative_pos_gene=="upstream" and gene_distance_class=="50kb"'),
        'up_lt_40kb': qcount('relative_pos_gene=="upstream" and gene_distance_class=="40kb"'),
        'up_lt_30kb': qcount('relative_pos_gene=="upstream" and gene_distance_class=="30kb"'),
        'up_lt_20kb': qcount('relative_pos_gene=="upstream" and gene_distance_class=="20kb"'),
        'up_lt_10kb': qcount('relative_pos_gene=="upstream" and gene_distance_class=="10kb"'),
        'in_Exon': qcount('In_gene_type == "Exon"'),
        'in_Intron': qcount('In_gene_type == "Intron"'),
        'down_lt_10kb': qcount('relative_pos_gene=="downstream" and gene_distance_class=="10kb"'),
        'down_lt_20kb': qcount('relative_pos_gene=="downstream" and gene_distance_class=="20kb"'),
        'down_lt_30kb': qcount('relative_pos_gene=="downstream" and gene_distance_class=="30kb"'),
        'down_lt_40kb': qcount('relative_pos_gene=="downstream" and gene_distance_class=="40kb"'),
        'down_lt_50kb': qcount('relative_pos_gene=="downstream" and gene_distance_class=="50kb"'),
        'down_gt_50kb': qcount('relative_pos_gene=="downstream" and gene_distance_class==">50kb"'),

        'in_oncogene': qcount('In_oncogene == "Yes"'),
        'in_suppressor': qcount('In_suppressor == "Yes"'),

        'up_lt_2.5kb': qcount('relative_pos_TSS=="upstream" and TSS_distance_class=="2.5kb"'),
        'up_lt_2.0kb': qcount('relative_pos_TSS=="upstream" and TSS_distance_class=="2.0kb"'),
        'up_lt_1.5kb': qcount('relative_pos_TSS=="upstream" and TSS_distance_class=="1.5kb"'),
        'up_lt_1.0kb': qcount('relative_pos_TSS=="upstream" and TSS_distance_class=="1.0kb"'),
        'up_lt_0.5kb': qcount('relative_pos_TSS=="upstream" and TSS_distance_class=="0.5kb"'),
        'in_TSS': qcount('relative_pos_TSS=="inner" and TSS_distance_class=="0.0kb"'),
        'down_lt_0.5kb': qcount('relative_pos_TSS=="downstream" and TSS_distance_class=="0.5kb"'),
        'down_lt_1.0kb': qcount('relative_pos_TSS=="downstream" and TSS_distance_class=="1.0kb"'),
        'down_lt_1.5kb': qcount('relative_pos_TSS=="downstream" and TSS_distance_class=="1.5kb"'),
        'down_lt_2.0kb': qcount('relative_pos_TSS=="downstream" and TSS_distance_class=="2.0kb"'),
        'down_lt_2.5kb': qcount('relative_pos_TSS=="downstream" and TSS_distance_class=="2.5kb"'),

        **{
            state: row_metric(
                mapping_site_feature.loc[mapping_site_feature['chromatin_state'] == state]
            )
            for state in chromatin_states
        },

        'far_away_TSS': qcount('far_away_TSS=="Yes"'),
        'far_away_cancerGene': qcount('far_away_cancerGene=="Yes"'),
        'far_away_miRNA': qcount('far_away_miRNA=="Yes"'),
        'out_gene': qcount('out_gene=="Yes"'),
        'out_conserved_region': qcount('out_conserved_region=="Yes"'),
        'combine_safe_harbor': qcount('combined_safe_harbor == True'),
    }

    total_events = (
        float(mapping_site_feature[weight_col].sum())
        if use_fragments
        else float(len(mapping_site_feature))
    )
    stat_mode = f'fragment_sum({weight_col})' if use_fragments else 'site_count(rows)'
    print(f'[cal_feature_ratio] file={site_with_feature_mark}\tmode={stat_mode}\ttotal_events={total_events:g}')
    print(feature_count)
    feature_count = pd.DataFrame([feature_count]).T 
    feature_ratio = (feature_count / total_events)

    feature_ratio = feature_ratio.reset_index()
    feature_ratio.columns = ['category', 'ratio']
    feature_count = feature_count.reset_index()
    feature_count.columns = ['category', 'count']
    return feature_ratio, feature_count, int(round(total_events)), stat_mode


def trans_zero_to_minor_positive(number):
    if number == 0:
        return 1 / 2
    else:
        return number


if __name__ == '__main__':
    params = parse_params()
    random = params.random
    outpath = params.outpath
    outlabel = params.outlabel
    outfile = os.path.join(outpath, outlabel + '.png')

    sample_file = {}
    for sample, file in zip(params.sample_label.split(','), params.mapping_file.split(',')):
        sample_file[sample] = file
    print(sample_file, random, outfile)

    try:
        with open(os.path.join(outpath, outlabel + '_event_number.yaml')) as fR:
            sample_event_number = yaml.safe_load(fR)
        Ratio = pd.read_csv(os.path.join(outpath, outlabel + '_ratio.csv'), sep='\t', index_col=0)
        Count = pd.read_csv(os.path.join(outpath, outlabel + '_count.csv'), sep='\t', index_col=0)
        odds_ratio = pd.read_csv(os.path.join(outpath, outlabel + '_odds_ratio.csv'), sep='\t', index_col=0)
        _smt = os.path.join(outpath, outlabel + '_fragment_stat_mode.tsv')
        if os.path.isfile(_smt):
            print('[Per-input counting mode (cached with ratio tables)]')
            print(pd.read_csv(_smt, sep='\t').to_string(index=False))

    except Exception as error:
        sample_event_number = {}
        stat_mode_rows = []
        backgroud_feature_ratio, backgroud_feature_count, bg_events, bg_stat_mode = cal_feature_ratio(random)
        stat_mode_rows.append(
            {'sample': '__BACKGROUD__', 'path': random, 'stat_mode': bg_stat_mode, 'total_events': bg_events}
        )
        backgroud_feature_ratio.columns = ['category', '__BACKGROUD__']
        backgroud_feature_count.columns = ['category', '__BACKGROUD__']
        mapping_feature_ratio = []
        mapping_feature_count = []
        for sample, mapping_file in sample_file.items():
            ratio, count, event_number, smode = cal_feature_ratio(mapping_file)
            if not (sample == 'backgroud' and mapping_file == random):
                stat_mode_rows.append(
                    {'sample': sample, 'path': mapping_file, 'stat_mode': smode, 'total_events': event_number}
                )
            sample_event_number[sample] = event_number
            mapping_feature_ratio.append(ratio)
            mapping_feature_ratio[-1].columns = ['category', sample]
            backgroud_feature_ratio = pd.merge(backgroud_feature_ratio, mapping_feature_ratio[-1], left_on=['category'], right_on=['category'], how='left')
            mapping_feature_count.append(count)
            mapping_feature_count[-1].columns = ['category', sample]
            backgroud_feature_count = pd.merge(backgroud_feature_count, mapping_feature_count[-1], left_on=['category'], right_on=['category'], how='left')

        with open(os.path.join(outpath, outlabel + '_event_number.yaml'), 'w') as fw:
            yaml.dump(sample_event_number, fw)
        stat_mode_tsv = os.path.join(outpath, outlabel + '_fragment_stat_mode.tsv')
        pd.DataFrame(stat_mode_rows).to_csv(stat_mode_tsv, sep='\t', index=False)
        print('[Per-input counting modes] written to', stat_mode_tsv)
        print(pd.DataFrame(stat_mode_rows).to_string(index=False))
        backgroud_feature_ratio.to_csv(os.path.join(outpath, outlabel + '_ratio.csv'), sep='\t')
        backgroud_feature_count.to_csv(os.path.join(outpath, outlabel + '_count.csv'), sep='\t') 

        Ratio = backgroud_feature_ratio.copy() 
        Count = backgroud_feature_count.copy() 
        for sample in sample_file:
            backgroud_feature_ratio[sample] = backgroud_feature_ratio[sample] / backgroud_feature_ratio['__BACKGROUD__']
        backgroud_feature_ratio['__BACKGROUD__'] = backgroud_feature_ratio['__BACKGROUD__']/backgroud_feature_ratio['__BACKGROUD__']

        backgroud_feature_ratio = backgroud_feature_ratio.set_index('category', drop=True)
        backgroud_feature_ratio = backgroud_feature_ratio.applymap(trans_zero_to_minor_positive)
        backgroud_feature_ratio = backgroud_feature_ratio.applymap(np.log2)
        backgroud_feature_ratio.to_csv(os.path.join(outpath, outlabel + '_odds_ratio.csv'), sep='\t')
        odds_ratio = backgroud_feature_ratio.copy() 

    content = []
    labels = [
        '>50', '50', '40', '30', '20', '10', 'In Exon', 'In Intron', '10', '20', '30', '40', '50', '>50',
        'In Oncogene',
        'In TSG',
        '2.5', '2.0', '1.5', '1.0', '0.5', 'TSS', '0.5', '1.0', '1.5', '2.0', '2.5',
        'Active TSS',
        'Flanking active TSS',
        "Transcr. at gene 5' and 3'",
        'Strong transcription',
        'Weak transcription',
        'Genic enhancers',
        'Enhancers',
        'ZNF genes + repeats',
        'Heterochromatin',
        'Bivalent/poised TSS',
        'Flanking bivalent TSS/Enh',
        'Bivalent enhancer',
        'Repressed Polycomb',
        'Weak repressed Polycomb',
        'Quiescent/low',
        '50kb from TSS',
        '300kb from cancer gene',
        '300kb from miRNA gene',
        'Outside of gene',
        'Outside of conserved region',
        'Combined',]

    xlabel_font = FontProperties( size=24)
    ylabel_font = FontProperties( size=24)
    title_font = FontProperties(size=26,weight='bold')

    clist = ['#006699', '#F8F8F8', '#d63600'] 
    newcmp = LinearSegmentedColormap.from_list('my_linear_color', clist)

    blocks = [[0, 14], [14, 16], [16, 27], [27, 42], [42, 48]]  
    row_num = len(sample_file)
    col_num = len(blocks)

    if os.path.isfile(params.sample_order):
        ylabels = pd.read_csv(params.sample_order,header=None)[0].to_list()
        ylabels = [i.replace('.', '-') for i in ylabels]
    else: 
        ylabels = list(sample_file.keys())
    print(ylabels)

    if len(ylabels) != row_num:
        info = 'ylabels is not matching row_num.'
        raise Exception(info)

    fig = plt.figure(figsize=(30, len(ylabels)), dpi=300)
    axes_colorbar_from = []

    colorbar_range = [[-4, 4]] * 5 
    _cbar = int(np.floor(odds_ratio.min().min()))
    cbar_ = int(np.ceil(odds_ratio.max().max()))
    colorbar_range = [[  _cbar ,  cbar_  ]] * 5  
    gs = gridspec.GridSpec(row_num, col_num, width_ratios=[b[1] - b[0] for b in blocks]) 
    titles = ['Gene', 'Cancer Gene', 'Transcript Start Site\n(TSS)', 'Chromatin state(ChromHMM)', 'Safe harbor']
    for i in range(row_num * col_num):

        axes = plt.subplot(gs[i])
        block = blocks[i % col_num]
        to_plot = odds_ratio[ylabels].T.iloc[[i // col_num], block[0]:block[1]]
        label = labels[block[0]:block[1]]
        sns.heatmap(to_plot,
                    square=False,
                    cmap=newcmp,
                    vmin = colorbar_range[i % col_num][0],
                    vmax = colorbar_range[i % col_num][1],
                    center = 0,
                    ax = axes,
                    cbar=False, linewidths=0.0, linecolor='white')
        axes_colorbar_from.append(axes) 
        if i // col_num == 0:  
            axes.set_xticks([])
            axes.set_xlabel('')
            axes.set_yticks([])
            axes.set_title(titles[i % col_num], fontproperties=title_font, y=1.2)

        elif i // col_num == row_num - 1: 
            axes.set_xticks([0.8 + j for j in range(to_plot.shape[1])])
            axes.set_xticklabels(label, rotation=70, ha='right', fontproperties=xlabel_font)
            axes.tick_params(axis='x', which='both', bottom=False, pad=20)

            if i % col_num in [0, 2]:
                axes.set_xlabel('(kb)', fontproperties=xlabel_font)
            else:
                axes.set_xlabel('')

            axes.set_yticks([])
        else:
            axes.set_xticks([])
            axes.set_xlabel('')
            axes.set_yticks([])

        if i % col_num == col_num - 1: 
            _n_disp = int(round(sample_event_number[ylabels[i // col_num]]))
            axes.set_ylabel(f'{ylabels[i // col_num]}\nn={_n_disp:,}', rotation=0, ha='left', va='center', fontstyle='italic', font=ylabel_font)
            axes.yaxis.set_label_coords(1.1, 0.5) 
    plt.subplots_adjust(wspace=0.04, hspace=0.2)

    cbar_ax = plt.axes((0.13, -3 / len(ylabels), 0.4, 0.5 / len(ylabels)))   
    cbar = fig.colorbar(axes_colorbar_from[0].collections[0], cax =cbar_ax, orientation ="horizontal")
    cbar.ax.tick_params(axis='x', labelsize=20, color='black', direction='in', pad=10)
    cbar.outline.set_edgecolor('white')
    cbar.ax.set_title('Enrichment Ratio', fontproperties=title_font, y=1.03,loc='left')
    x_ticks = [i for i in range(_cbar,cbar_+1)]
    cbar.ax.set_xticks(x_ticks)
    xticklabels = []
    for i in x_ticks:
        if i < 0:
            xticklabels.append(f'1/{2**abs(i)}')
        elif i == 0:
            xticklabels.append(f'1')
        else:
            xticklabels.append(f'{2**i}')
    cbar.ax.set_xticklabels(xticklabels)

    plt.savefig(outfile, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(outfile+'.pdf', dpi=300, bbox_inches='tight', transparent=True)