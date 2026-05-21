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
#   Get a representative peptide length from a set of sequences.                            #
# ------------------------------------------------------------------------------------------#

import numpy as np
import pandas as pd


def get_length_common( sequences, step = 4 ):
    sequences = np.array( sequences )
    temp_size = len( sequences )
    temp_seq_len   = [ len(sequences[i]) for i in range(temp_size) ]
    temp_seq_index = [ i for i in range(temp_size) ]
    temp_seq_data  = pd.DataFrame( {"seq_index":temp_seq_index, "seq_len":temp_seq_len } ).sort_values( ['seq_len'], ascending = True )
    temp_seq_data_1  = temp_seq_data.copy()
    temp_seq_data_2  = temp_seq_data.copy()
    temp_seq_data_1['len_group']  = temp_seq_data_1.apply( lambda x:"m" + str( (x['seq_len'] + 0) // 4 ) , axis = 1)
    temp_seq_data_2['len_group']  = temp_seq_data_2.apply( lambda x:"n" + str( (x['seq_len'] + 2) // 4 ) , axis = 1)
    temp_seq_data  = pd.concat( [temp_seq_data_1, temp_seq_data_2] )
    temp_seq_data_grouped = temp_seq_data.groupby( ['len_group'] ).size().reset_index( name = "size" ).merge( temp_seq_data.groupby( ['len_group'] )['seq_len'].median().reset_index(), left_on="len_group", right_on="len_group" ).sort_values( ['size', 'seq_len'], ascending = [False, False] )
    temp_seq_data_grouped_top  = temp_seq_data_grouped['len_group'].iloc[0]
    temp_seq_data_group_choose = temp_seq_data.query( 'len_group == @temp_seq_data_grouped_top' )
    pass
    temp_seq_data_group_choose_top = np.bincount( temp_seq_data_group_choose['seq_len'] ).argmax()
    temp_seq_data_group_choose_top_index = temp_seq_data_group_choose.query( 'seq_len == @temp_seq_data_group_choose_top' )['seq_index'].iloc[0]
    return( temp_seq_data_group_choose_top_index, sequences[ temp_seq_data_group_choose_top_index ] )
