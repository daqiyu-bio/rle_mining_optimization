#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------------------#
# Copyright (c) 2026 Beijing Astragenomics Technology Co., Ltd. All rights reserved.        #
# Licensing: see the LICENSE file in the repository root.                                   #
# version : 1.0                                                                             #
#                                                                                           #
# Author:                                                                                   #
#   Daqi Yu <daqi.yu@astragenomics.com>                                                     #
#                                                                                           #
# Description:                                                                              #
#   Filter and reconstruct annotation records for non-LTR loci.                             #
# ------------------------------------------------------------------------------------------#

import re
import numpy as np
import pandas as pd


def get_orf_pfam_nonLTR_locus( item_fna_orf_pfam ):
    item_fna_orf_pfam = item_fna_orf_pfam.query( 'bit_score > 10' )
    item_fna_orf_pfam_RVT1  = item_fna_orf_pfam[ item_fna_orf_pfam['domain_name'].isin( ['RVT_1'] )   ][ ['genome_ID', 'chr', 'orf_strand', 'orf_start', 'orf_end', 'orf_ID'] ].drop_duplicates().copy()
    item_fna_orf_pfam_RLE   = item_fna_orf_pfam[ item_fna_orf_pfam['domain_name'].str.endswith('RLE') ][ ['genome_ID', 'chr', 'orf_strand', 'orf_start', 'orf_end', 'orf_ID'] ].drop_duplicates().copy()
    item_fna_orf_pfam_RVT1['orf_type']  = 'nonLTR'
    if item_fna_orf_pfam_RVT1.shape[0] > 0:
        pass
    else:
        return( None )
    item_fna_orf_pfam_parse_locus = item_fna_orf_pfam_RVT1.sort_values( ['genome_ID', 'chr', 'orf_start'] ).copy()
    item_fna_orf_pfam_parse_locus['locus_type'] = 'nonLTR'
    item_fna_orf_pfam_parse_locus['locus_ID']   = [ i for i in range(item_fna_orf_pfam_parse_locus.shape[0]) ]
    item_fna_orf_pfam_parse_locus_res = item_fna_orf_pfam[ ['genome_ID', 'chr', 'orf_start', 'orf_end', 'orf_strand', 'orf_ID', 'orf_seq', 'orf_pep', 'orf_pep_len'] ].drop_duplicates().merge( item_fna_orf_pfam_parse_locus[ ['genome_ID', 'chr', 'orf_ID', 'orf_type', 'locus_type', 'locus_ID'] ], left_on = ['genome_ID', 'chr', 'orf_ID'], right_on = ['genome_ID', 'chr', 'orf_ID'] ).sort_values( ['locus_ID', 'orf_type'] ).copy()
    item_fna_orf_pfam_parse_locus_res.loc[item_fna_orf_pfam_parse_locus_res['locus_type'].isin( ['nonLTR'] ), 'locus_group'] = "non_LTR"
    item_fna_orf_pfam_parse_locus_res['locus_group'] = item_fna_orf_pfam_parse_locus_res.apply(lambda x:x['locus_group'] if x['locus_group'] == x['locus_group'] else "others", axis=1)
    if item_fna_orf_pfam_parse_locus_res.shape[0] == 0:
        return( None )
    item_fna_orf_pfam_parse_locus_res.loc[item_fna_orf_pfam_parse_locus_res['orf_ID'].isin( item_fna_orf_pfam_RVT1['orf_ID'] ), 'RVT_1'] = "Y"
    item_fna_orf_pfam_parse_locus_res.loc[item_fna_orf_pfam_parse_locus_res['orf_ID'].isin( item_fna_orf_pfam_RLE['orf_ID'] ), 'RLE']   = "Y"
    item_fna_orf_pfam_parse_locus_res['locus_strand']  = item_fna_orf_pfam_parse_locus_res.apply(lambda x:x['orf_strand'], axis=1)
    item_fna_orf_pfam_parse_locus_res['locus_left_end']  = item_fna_orf_pfam_parse_locus_res.apply(lambda x:x['orf_start'], axis=1)
    item_fna_orf_pfam_parse_locus_res['locus_right_end'] = item_fna_orf_pfam_parse_locus_res.apply(lambda x:x['orf_end'], axis=1)
    return( item_fna_orf_pfam_parse_locus_res.drop(["locus_ID"], axis=1) )
