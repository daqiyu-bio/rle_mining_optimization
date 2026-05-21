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
#   MSA-based consensus parsing for RLE/non-LTR workflows.                                  #
# ------------------------------------------------------------------------------------------#

import re
import numpy as np
from Bio import Align
from Bio import AlignIO


# Recommended:
#aligner_strict = Align.PairwiseAligner()
#aligner_strict.mode = 'local'
#aligner_strict.match_score = 2
#aligner_strict.mismatch_score = -12
#aligner_strict.open_gap_score = -12
#aligner_strict.extend_gap_score = -1

def align_consensus_parse_end( msa_alignment, aligner_strict, align_final_choice, conserve = True, mismatch_size = 40, msa_end = "right"  ):
    align_size     = len( msa_alignment )
    align_len      = len( msa_alignment[0,:] )
    """
    ###Step1. only some sequences were chosen to be kept for consensus sequence construction.
    """    
    assert align_size == len( align_final_choice )
    align_size_final  = sum( align_final_choice )
    if align_size_final <= 3:
        identity = 2/3
    elif conserve == True:
        identity = 0.75
    else:
        identity = 0.66
    pass
    align_consensus  = ""
    for i in range( align_len ):
        site_freq = np.unique( np.array( list( msa_alignment[:,i] ) )[align_final_choice], return_counts = True )
        site_max  = np.max( site_freq[1] )
        if site_max/align_size_final >= identity:
            site_base = site_freq[0][ np.argmax( site_freq[1] ) ]
        else:
            i_random  = i % align_size_final
            site_base = np.array( list( msa_alignment[:,i] ) )[align_final_choice][ i_random ]
        align_consensus = align_consensus + site_base
    pass
    align_consensus_start = []
    align_consensus_end   = []
    for i in range( align_size ):
        if align_final_choice[i] == True:
            align_res = aligner_strict.align( re.sub( "-", "", align_consensus ), re.sub( "-", "", str(msa_alignment[i].seq) ) )
            align_consensus_start.append( align_res[0].coordinates.T[0][0] )
            align_consensus_end.append( align_res[0].coordinates.T[-1][0] )
    """
    End choice difference exits when dealing with 5' split MSA and 3' split MSA.
    """
    align_consensus_start_ref    = np.median( align_consensus_start )
    align_consensus_end_ref      = np.median( align_consensus_end )
    if msa_end == "right":
        align_consensus_choice_start = list( np.array(align_consensus_start ) - align_consensus_start_ref  <= mismatch_size )
        align_consensus_choice_end   = [ True for i in range(align_size_final) ]
    else:
        align_consensus_choice_start = [ True for i in range(align_size_final) ]
        align_consensus_choice_end   = list( align_consensus_end_ref  -  np.array(align_consensus_end  )   <= mismatch_size )
    align_consensus_start        = np.array( align_consensus_start )[ align_consensus_choice_start ]
    align_consensus_end          = np.array( align_consensus_end )[ align_consensus_choice_end ]
    align_consensus_choice_start = np.array( [ align_consensus_choice_start.pop(0) if align_final_choice[i] == True else align_final_choice[i] for i in range(align_size)  ] )
    align_consensus_choice_end   = np.array( [ align_consensus_choice_end.pop(0)   if align_final_choice[i] == True else align_final_choice[i] for i in range(align_size)  ] )
    align_consensus_choice       = align_final_choice & align_consensus_choice_start & align_consensus_choice_end
    pass
    align_size_consensus      = sum( align_consensus_choice )
    """
    ####Setp3 In some cases, we filtered too many aligned sequence, we then tried to re-built the consensus region.
    """
    if align_size_consensus/align_size_final < identity:
        align_consensus  = ""
        for i in range( align_len ):
            site_freq = np.unique( np.array( list( msa_alignment[:,i] ) )[align_consensus_choice], return_counts = True )
            site_max  = np.max( site_freq[1] )
            if site_max/align_size_consensus >= identity:
                site_base = site_freq[0][ np.argmax( site_freq[1] ) ]
            else:
                i_random  = i % align_size_consensus
                site_base = np.array( list( msa_alignment[:,i] ) )[align_consensus_choice][ i_random ]
            align_consensus = align_consensus + site_base
        align_consensus_start = []
        align_consensus_end   = []
        for i in range( align_size ):
            if align_consensus_choice[i] == True:
                align_res = aligner_strict.align( re.sub( "-", "", align_consensus ), re.sub( "-", "", str(msa_alignment[i].seq) ) )
                align_consensus_start.append( align_res[0].coordinates.T[0][0] )
                align_consensus_end.append( align_res[0].coordinates.T[-1][0] )
    pass
    align_consensus_pos = []
    for i in range( align_len ):
        if align_consensus[i] != "-":
            align_consensus_pos.append( i )
    align_consensus_msa_start = align_consensus_pos[ int( np.ceil( np.median( align_consensus_start ) ) ) ]
    align_consensus_msa_end   = align_consensus_pos[ int( np.floor( np.median( align_consensus_end ) ) ) - 1 ]
    pass
    align_consensus = ""
    for i in range( align_consensus_msa_start, align_consensus_msa_end + 1 ):
        site_freq = np.unique( np.array( list( msa_alignment[:, i] ) )[align_consensus_choice], return_counts = True )
        site_base = site_freq[0][ np.argmax( site_freq[1] ) ]
        align_consensus = align_consensus + site_base
    pass
    align_consensus_choice_seq   = re.sub("-", "", align_consensus).upper()
    return( align_consensus_choice, align_consensus_choice_seq, align_consensus_msa_start, align_consensus_msa_end  )


def align_pairwise_coverage( msa_alignment, i, j, msa_end = "right" , match_score = 2, mismatch_score = -10, open_gap_score = -10, extend_gap_score = -1):
    align_len   = len( msa_alignment[0] )
    seq_len_1   = len( str(msa_alignment[i].seq).replace( '-', '') )
    seq_len_2   = len( str(msa_alignment[j].seq).replace( '-', '') )
    align_score = []    
    k_previous  = ''
    align_range = [i for i in range(align_len)]
    if msa_end == "right":
        align_range = align_range[::-1]
    else:
        pass
    match_count = 0
    for k in align_range:
        if msa_alignment[i][k] == msa_alignment[j][k]:
            if msa_alignment[i][k] == '-':
                k_score = -0.001
            else:
                k_score = match_score
                k_previous = "non_gap"
                match_count += 1
        else:
            if msa_alignment[i][k] == '-' or msa_alignment[j][k] == '-':
                if k_previous == "gap":
                    k_score = extend_gap_score
                else:
                    k_score = open_gap_score
                k_previous = "gap"
            else:
                k_score = mismatch_score
                k_previous = 'non_gap'
        align_score.append( k_score )
        if match_count == 10:
            align_score = [ 0 for i in range( len(align_score) ) ]
    align_cum_score = np.cumsum( align_score )
    align_cum_max   = np.argmax( align_cum_score )
    if msa_end == "right":
        align_span_1     = len( str(msa_alignment[i].reverse_complement()[0:align_cum_max].seq).replace( '-', '') ) + 1
        align_span_2     = len( str(msa_alignment[j].reverse_complement()[0:align_cum_max].seq).replace( '-', '') ) + 1
    else:
        align_span_1     = len( str(msa_alignment[i][0:align_cum_max].seq).replace( '-', '') ) + 1
        align_span_2     = len( str(msa_alignment[j][0:align_cum_max].seq).replace( '-', '') ) + 1
    if seq_len_1 == 0 or seq_len_2 == 0:
        return( 0, 0, seq_len_1, seq_len_2, align_len )
    else:
        return( align_span_1/seq_len_1, align_span_2/seq_len_2, seq_len_1 - align_span_1, seq_len_2 - align_span_2, align_len - align_cum_max )


def align_truncation_boundary_pos( msa_alignment, match_score = 2, mismatch_score = -10, open_gap_score = -10, extend_gap_score = -1 ):
    align_size     = len( msa_alignment )
    align_len      = len( msa_alignment[0,:] )
    align_identical = []
    for i in range( align_len ):
        site_freq = np.unique( np.array( list( msa_alignment[:,i] ) ), return_counts = True )
        site_max  = np.max( site_freq[1] )
        site_base = site_freq[0][ np.argmax( site_freq[1] ) ]
        if site_base != "-" and site_max/align_size >= 0.8:
            align_identical.append( i )
    align_identical_base = int( np.median( align_identical ) )
    msa_left  = msa_alignment[:, 0:align_identical_base]
    msa_right = msa_alignment[:, align_identical_base:]
    msa_raw_left  = [ str(msa_left[i].seq).replace('-', '').upper() for i in range(align_size) ]
    msa_raw_right = [ str(msa_right[i].seq).replace('-', '').upper() for i in range(align_size) ]
    msa_raw_left_len  = [ len(msa_raw_left[i]) for i in range(align_size) ] 
    msa_raw_right_len = [ len(msa_raw_right[i]) for i in range(align_size) ]
    """
    Step1. Mark duplicate sequence with Kmer. That is, 12-mer redundancy at first or last 200 bp. Ns are excluded.
    """
    msa_dedup_left    = [ True ] *  align_size
    msa_dedup_right   = [ True ] *  align_size
    msa_all_Kmer   = set( [] )
    for i in np.argsort( msa_raw_left_len )[::-1]:
        align_seq_Kmer = [ msa_raw_left[i][k:(k + 12)] for k in range( len(msa_raw_left[i]) - 11 ) if msa_raw_left[i][k] != "N" and msa_raw_left[i][k+11] != "N" ]
        align_seq_left_Kmer   = set( align_seq_Kmer[0:200] )
        if len(align_seq_left_Kmer) > 0 and len( align_seq_left_Kmer.intersection( msa_all_Kmer) )/len(align_seq_left_Kmer) < 0.4:
            msa_all_Kmer  = msa_all_Kmer.union( set( align_seq_Kmer ) )
            msa_dedup_left[i] = True
        else:
            msa_dedup_left[i] = False
    msa_all_Kmer   = set( [] )
    for i in np.argsort( msa_raw_right_len )[::-1]:
        align_seq_Kmer = [ msa_raw_right[i][k:(k + 12)] for k in range( len(msa_raw_right[i]) - 11 ) if msa_raw_right[i][k] != "N" and msa_raw_right[i][k+11] != "N" ]
        align_seq_right_Kmer  = set( align_seq_Kmer[-200:] )
        if len(align_seq_right_Kmer) > 0 and len( align_seq_right_Kmer.intersection( msa_all_Kmer) )/len(align_seq_right_Kmer) < 0.4:
            msa_all_Kmer  = msa_all_Kmer.union( set( align_seq_Kmer ) )
            msa_dedup_right[i] = True
        else:
            msa_dedup_right[i] = False
    pass
    """
    Step2. Mark 5' truncated sequence.
    """
    align_left_pos       = [ dict() for i in range(align_size) ]
    align_right_pos      = [ dict() for i in range(align_size) ]
    align_left_msa_pos   = [ dict() for i in range(align_size) ]
    align_right_msa_pos  = [ dict() for i in range(align_size) ]
    for i in range( align_size ):
        for j in range( i + 1, align_size):
            align_span = align_pairwise_coverage( msa_left, i, j, msa_end = "right",  match_score = match_score, mismatch_score = mismatch_score, open_gap_score = open_gap_score, extend_gap_score = extend_gap_score )
            align_left_pos[i][j]     = align_span[2]
            align_left_pos[j][i]     = align_span[3]
            align_left_msa_pos[i][j] = align_span[4]
            align_left_msa_pos[j][i] = align_span[4]
    for i in range( align_size ):
        for j in range( i + 1, align_size):
            align_span = align_pairwise_coverage( msa_right, i, j, msa_end = "left",  match_score = match_score, mismatch_score = mismatch_score, open_gap_score = open_gap_score, extend_gap_score = extend_gap_score )
            align_right_pos[i][j]     = align_span[2]
            align_right_pos[j][i]     = align_span[3]
            align_right_msa_pos[i][j] = align_span[4]
            align_right_msa_pos[j][i] = align_span[4]
    pass
    msa_dedup_left_size = sum( msa_dedup_left )
    """
    At least two unique sequences are required to define 5' truncation.
    """
    if msa_dedup_left_size > 1:
        msa_left_msa_pos  = [ min([align_left_msa_pos[i][j] for j in align_left_msa_pos[i].keys() if msa_dedup_left[j] == True ]) for i in range( align_size ) ]
        msa_pos_5         = np.median( [ msa_left_msa_pos[i] for i in range( align_size ) if msa_dedup_left[i] == True ] )
        msa_5_truncation  = np.array( msa_left_msa_pos - msa_pos_5 >= 50 )
    else:
        msa_5_truncation  = np.array( [ False ] * align_size )
    return( align_identical_base, msa_dedup_left, msa_dedup_right, msa_5_truncation )




def align_consensus_parse_RNA( msa_path, aligner_strict, conserve = True, mismatch_size = 40 ):
    msa_alignment  = AlignIO.read(open( msa_path ), "fasta")
    align_size     = len( msa_alignment )
    align_len      = len( msa_alignment[0,:] )
    """
    ###Step1. descarding the potential 5' truncated locus and duplicated locus
    """    
    align_identical_base, msa_dedup_left, msa_dedup_right, msa_5_truncation  = align_truncation_boundary_pos( msa_alignment, match_score = 2, mismatch_score = -10, open_gap_score = -10, extend_gap_score = -1 )
    msa_not_truncated        = np.array( msa_5_truncation != True )
    align_final_choice_left  = np.array( msa_dedup_left ) & np.array( msa_not_truncated )
    align_final_choice_right = np.array( msa_dedup_right )
    """
    ###Step2. Generating 5' and 3' consensus sequence with the chosen sequences
    """    
    align_identical = []
    for i in range( align_len ):
        site_freq = np.unique( np.array( list( msa_alignment[:,i] ) ), return_counts = True )
        site_max  = np.max( site_freq[1] )
        site_base = site_freq[0][ np.argmax( site_freq[1] ) ]
        if site_base != "-" and site_max/align_size >= 0.8:
            align_identical.append( i )
    align_identical_base = int( np.median( align_identical ) )  ###Using the newly generated split postition
    msa_left  = msa_alignment[:, 0:align_identical_base]
    msa_right = msa_alignment[:, align_identical_base:]
    res1   = align_consensus_parse_end( msa_left,  aligner_strict = aligner_strict, align_final_choice = align_final_choice_left, conserve = conserve, mismatch_size = mismatch_size, msa_end = "right" )
    res2   = align_consensus_parse_end( msa_right, aligner_strict = aligner_strict, align_final_choice = align_final_choice_right, conserve = conserve, mismatch_size = mismatch_size, msa_end = "left" )
    if res1 is None or res2 is None:
        return( None )
    else:
        align_consensus_choice     = (res1[0], res2[0])
        align_consensus_choice_seq = res1[1] + res2[1]
        align_consensus_msa_start  = res1[2]
        align_consensus_msa_end    = align_identical_base + res2[3]
        return( align_consensus_choice, align_consensus_choice_seq, align_consensus_msa_start, align_consensus_msa_end )
