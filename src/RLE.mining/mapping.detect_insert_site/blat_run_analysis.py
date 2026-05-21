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
#   Post-process BLAT PSL (block / indel filters).                                          #
# ------------------------------------------------------------------------------------------#

import os
import re
import pandas as pd
from Bio import SeqIO


def blat_block_parse( query_block_start, target_block_start, block_size, max_gap_span = 200 ):
    """
    Note that in BLAT blocks, reverse-complement coordinates are used on the query strand.
    The target coordinates are always on the plus strand.
    """
    block_query_span    = [(int(i), int(i) + int(k)) for i,k in zip(query_block_start,  block_size) if re.match("[0-9]", i) and re.match("[0-9]", k)]
    block_target_span   = [(int(j), int(j) + int(k)) for j,k in zip(target_block_start, block_size) if re.match("[0-9]", j) and re.match("[0-9]", k)]
    block_query_update  = []
    block_target_update = []
    item_query_start    = block_query_span[0][0]
    item_query_end      = block_query_span[0][1]
    item_target_start   = block_target_span[0][0]
    item_target_end     = block_target_span[0][1]
    block_count_total   = len(block_target_span)
    if block_count_total == 0:
        return( None )
    elif block_count_total == 1:
        block_query_update  = block_query_span
        block_target_update = block_target_span
    else:
        for i in range( 0, block_count_total - 1 ):
            if block_target_span[i+1][0] - item_target_end < max_gap_span and block_query_span[i+1][0] - item_query_end < max_gap_span:
                item_query_end    = block_query_span[i+1][1]
                item_target_end   = block_target_span[i+1][1]
            else:
                block_query_update.append(  (item_query_start, item_query_end) )
                block_target_update.append( (item_target_start, item_target_end) )
                item_query_start  = block_query_span[i+1][0]
                item_query_end    = block_query_span[i+1][1]
                item_target_start = block_target_span[i+1][0]
                item_target_end   = block_target_span[i+1][1]
        block_query_update.append(  (item_query_start, item_query_end) )
        block_target_update.append( (item_target_start, item_target_end) )
    block_updated_count     = len( block_target_update )
    block_query_start_list  = [ block_query_update[i][0] for i in range(block_updated_count) ]
    block_query_end_list    = [ block_query_update[i][1] for i in range(block_updated_count) ]
    block_target_start_list = [ block_target_update[i][0] for i in range(block_updated_count) ]
    block_target_end_list   = [ block_target_update[i][1] for i in range(block_updated_count) ]
    block_updated_DataFrame =  pd.DataFrame( [block_query_start_list, block_query_end_list, block_target_start_list, block_target_end_list] ).T
    return( block_updated_DataFrame )


def blat_match_block_parse( blat_run_file ):
    blat_out  = pd.read_csv( blat_run_file, index_col=False, skiprows=5, sep='\s+', names=['match', 'mismatch', 'rep_match', 'N', 'query_gap_count', 'query_gap_base' ,'target_gap_count' ,'target_gap_base', 'strand', 'query',  'query_size', 'query_start', 'query_end', 'target',  'target_size', 'target_start', 'target_end', 'block_count', 'block_size', "block_query_start", "block_target_start" ])
    if blat_out.shape[0] > 0:
        pass
    else:
        print( "Blat result is empty." )
        return( None )
    blat_out_parse     = blat_out.copy()
    blat_out_parse['identity'] = blat_out_parse.apply( lambda x: float( '{0:.3f}'.format( x['match']/(x['match'] + x['mismatch']) ) ), axis=1)
    blat_out_DataFrame = pd.DataFrame()
    blat_out_DataFrame_columns = ['query_start', 'query_end', 'target_start', 'target_end', 'query', 'query_size', 'target', 'strand', 'identity','source_label']
    for i in range( blat_out_parse.shape[0] ):
        query_block_start  = blat_out_parse.iloc[i]['block_query_start'].split(',')
        target_block_start = blat_out_parse.iloc[i]['block_target_start'].split(',')
        block_size         = blat_out_parse.iloc[i]['block_size'].split(',')
        block_DataFrame    = blat_block_parse( query_block_start, target_block_start, block_size )
        block_DataFrame['query']  = blat_out_parse.iloc[i]['query']
        block_DataFrame['query_size']  = blat_out_parse.iloc[i]['query_size']
        block_DataFrame['target'] = blat_out_parse.iloc[i]['target']
        block_DataFrame['strand'] = blat_out_parse.iloc[i]['strand']
        block_DataFrame['identity'] = blat_out_parse.iloc[i]['identity']
        block_DataFrame['source_label'] = '{chrom}:{start}-{end}'.format(chrom=str(blat_out_parse.iloc[i]['target']),start= blat_out_parse.iloc[i]['target_start'],end = blat_out_parse.iloc[i]['target_end'])  # genomic locus label for this hit
        blat_out_DataFrame = pd.concat( [ blat_out_DataFrame, block_DataFrame ] )
    blat_out_DataFrame.columns = blat_out_DataFrame_columns
    blat_out_DataFrame = blat_out_DataFrame.sort_values(['query', 'strand', 'query_start'])
    return( blat_out_DataFrame )


def blat_match_block_indels( reference_fna_path, genome_fna_dir, blat_run_dir,  run_id = "temp", min_match_size = 20, min_del_size = 20, min_gap_size = 2000, max_gap_size = 7000, min_overlap_size = 200 ):
    """
    Please note the 0-based nature of blat output
    """
    blat_run_file = os.path.join( blat_run_dir, run_id + ".psl" )
    blat_out_DataFrame = blat_match_block_parse( blat_run_file )
    if blat_out_DataFrame is None:
        return( pd.DataFrame(), pd.DataFrame(), pd.DataFrame() )
    blat_out_DataFrame = blat_out_DataFrame[blat_out_DataFrame.apply( lambda x:x['target_end'] - x['target_start'], axis=1) > min_match_size]
    blat_out_DataFrame = blat_out_DataFrame.sort_values(['query', 'target', 'strand', 'query_start'], ascending = [True, True, True, True])
    blat_out_DataFrame_count = blat_out_DataFrame.shape[0]
    blat_out_DataFrame['ID'] = [i for i in range(blat_out_DataFrame_count)]
    pass
    blat_out_DataFrame_gaps = pd.DataFrame()
    blat_out_DataFrame_dels = pd.DataFrame()
    blat_out_DataFrame_gaps_columns  = ['query', 'query_size', 'target', 'strand', 'query_start_5', 'query_end_5', 'target_start_5', 'target_end_5','source_label_5', 'query_start_3', 'query_end_3', 'target_start_3', 'target_end_3' ,'source_label_3',]
    blat_out_DataFrame_dels_columns  = ['query', 'query_size', 'target', 'strand', 'query_start_5', 'query_end_5', 'target_start_5', 'target_end_5','source_label_5', 'query_start_3', 'query_end_3', 'target_start_3', 'target_end_3' ,'source_label_3',]
    if blat_out_DataFrame.shape[0] > 1:
        blat_out_DataFrame_pair   = blat_out_DataFrame.merge( blat_out_DataFrame, left_on = ['query', 'target', 'strand'], right_on = ['query', 'target', 'strand'] ).query( 'ID_x != ID_y' )
        blat_out_DataFrame_nested = blat_out_DataFrame_pair.query( 'query_start_x >= query_start_y and query_end_x <= query_end_y and target_start_x >= target_start_y and target_end_x <= target_end_y ' ).copy()
        blat_out_DataFrame_pair   = blat_out_DataFrame_pair[ - blat_out_DataFrame_pair['ID_x'].isin( blat_out_DataFrame_nested['ID_x'] ) ]
        blat_out_DataFrame        = blat_out_DataFrame[ - blat_out_DataFrame['ID'].isin( blat_out_DataFrame_nested['ID_x'] ) ]
    else:
        blat_out_DataFrame_pair   = pd.DataFrame()
    if blat_out_DataFrame_pair.shape[0] > 0:
        blat_out_DataFrame_gaps = blat_out_DataFrame_pair.query( 'query_start_y > query_start_x and query_end_y > query_end_x and query_start_y - query_end_x > 0 - @min_overlap_size and query_start_y - query_end_x < @min_overlap_size and target_start_y - target_end_x > @min_gap_size and target_start_y - target_end_x < @max_gap_size' )
        blat_out_DataFrame_gaps = blat_out_DataFrame_gaps[ ['query', 'query_size_x', 'target', 'strand', 'query_start_x', 'query_end_x', 'target_start_x', 'target_end_x','source_label_x', 'query_start_y', 'query_end_y', 'target_start_y', 'target_end_y' ,'source_label_y'] ]
        blat_out_DataFrame_dels = blat_out_DataFrame_pair.query( 'query_start_y > query_start_x and query_end_y > query_end_x and target_start_y - target_end_x > 0 - @min_overlap_size and target_start_y - target_end_x < @min_overlap_size and query_start_y - query_end_x - (target_start_y - target_end_x) > @min_del_size and query_start_y - query_end_x < query_size_x' )
        blat_out_DataFrame_dels = blat_out_DataFrame_dels[ ['query', 'query_size_x', 'target', 'strand', 'query_start_x', 'query_end_x', 'target_start_x', 'target_end_x','source_label_x', 'query_start_y', 'query_end_y', 'target_start_y', 'target_end_y','source_label_y' ] ]

    pass
    with open( reference_fna_path, 'rt') as fin:
        query_fna  = dict( [(sequence.id, sequence) for sequence in SeqIO.parse(fin, 'fasta')] )
    with open( genome_fna_dir, 'rt') as fin:
        genome_fna = dict( [(sequence.id, sequence) for sequence in SeqIO.parse(fin, 'fasta')] )
    if blat_out_DataFrame_gaps.shape[0] > 0:
        blat_out_DataFrame_gaps.columns  = blat_out_DataFrame_gaps_columns
        blat_out_DataFrame_gaps_f = blat_out_DataFrame_gaps.query('strand == "+"' ).copy()
        blat_out_DataFrame_gaps_b = blat_out_DataFrame_gaps.query('strand == "-"' ).copy()
        if blat_out_DataFrame_gaps_f.shape[0] > 0:
            blat_out_DataFrame_gaps_f['target_sequence']  =  blat_out_DataFrame_gaps_f.apply( lambda x: str(genome_fna[ x['target'] ][(x['target_start_5']):x['target_end_3']].seq ).upper(), axis=1 )
            blat_out_DataFrame_gaps_f['insert_sequence']  =  blat_out_DataFrame_gaps_f.apply( lambda x: str(genome_fna[ x['target'] ][(x['target_end_5']):x['target_start_3']].seq ).upper(), axis=1 )
        if blat_out_DataFrame_gaps_b.shape[0] > 0:
            blat_out_DataFrame_gaps_b['target_sequence']  =  blat_out_DataFrame_gaps_b.apply( lambda x: str(genome_fna[ x['target'] ][(x['target_start_5']):x['target_end_3']].reverse_complement().seq ).upper(), axis=1 )
            blat_out_DataFrame_gaps_b['insert_sequence']  =  blat_out_DataFrame_gaps_b.apply( lambda x: str(genome_fna[ x['target'] ][(x['target_end_5']):x['target_start_3']].reverse_complement().seq ).upper(), axis=1 )
            blat_out_DataFrame_gaps_b = blat_out_DataFrame_gaps_b.assign( new_query_start_5 = lambda x:x['query_size'] - x['query_end_5'], new_query_end_5 = lambda x:x['query_size'] - x['query_start_5'], new_query_start_3 = lambda x:x['query_size'] - x['query_end_3'], new_query_end_3 = lambda x:x['query_size'] - x['query_start_3']  )
            blat_out_DataFrame_gaps_b = blat_out_DataFrame_gaps_b.assign(  query_start_5 = lambda x:x['new_query_start_5'], query_end_5 = lambda x:x['new_query_end_5'], query_start_3 = lambda x:x['new_query_start_3'], query_end_3 = lambda x:x['new_query_end_3'] )
            blat_out_DataFrame_gaps_b = blat_out_DataFrame_gaps_b.drop( ["new_query_start_5", "new_query_end_5", "new_query_start_3", "new_query_end_3"] , axis=1)
        blat_out_DataFrame_gaps = pd.concat( [ blat_out_DataFrame_gaps_f, blat_out_DataFrame_gaps_b ] )
    pass
    if blat_out_DataFrame_dels.shape[0] > 0:
        blat_out_DataFrame_dels.columns  = blat_out_DataFrame_dels_columns
        blat_out_DataFrame_dels_f = blat_out_DataFrame_dels.query('strand == "+"' ).copy()
        blat_out_DataFrame_dels_b = blat_out_DataFrame_dels.query('strand == "-"' ).copy()
        if blat_out_DataFrame_dels_f.shape[0] > 0:
            blat_out_DataFrame_dels_f['target_sequence']  =  blat_out_DataFrame_dels_f.apply( lambda x: str(genome_fna[ x['target'] ][(x['target_start_5']):x['target_end_3']].seq ).upper(), axis=1 )
            blat_out_DataFrame_dels_f['delete_sequence']  =  blat_out_DataFrame_dels_f.apply( lambda x: str(query_fna[ x['query'] ][(x['query_end_5']):x['query_start_3']].seq ).upper(), axis=1 )
        if blat_out_DataFrame_dels_b.shape[0] > 0:
            blat_out_DataFrame_dels_b['target_sequence']  =  blat_out_DataFrame_dels_b.apply( lambda x: str(genome_fna[ x['target'] ][(x['target_start_5']):x['target_end_3']].reverse_complement().seq ).upper(), axis=1 )
            blat_out_DataFrame_dels_b['delete_sequence']  =  blat_out_DataFrame_dels_b.apply( lambda x: str(query_fna[ x['query'] ].reverse_complement()[(x['query_end_5']):x['query_start_3']].reverse_complement().seq ).upper(), axis=1 )
            blat_out_DataFrame_dels_b = blat_out_DataFrame_dels_b.assign( new_query_start_5 = lambda x:x['query_size'] - x['query_end_5'], new_query_end_5 = lambda x:x['query_size'] - x['query_start_5'], new_query_start_3 = lambda x:x['query_size'] - x['query_end_3'], new_query_end_3 = lambda x:x['query_size'] - x['query_start_3']  )
            blat_out_DataFrame_dels_b = blat_out_DataFrame_dels_b.assign(  query_start_5 = lambda x:x['new_query_start_5'], query_end_5 = lambda x:x['new_query_end_5'], query_start_3 = lambda x:x['new_query_start_3'], query_end_3 = lambda x:x['new_query_end_3'] )
            blat_out_DataFrame_dels_b = blat_out_DataFrame_dels_b.drop( ["new_query_start_5", "new_query_end_5", "new_query_start_3", "new_query_end_3"] , axis=1)
        blat_out_DataFrame_dels = pd.concat( [ blat_out_DataFrame_dels_f, blat_out_DataFrame_dels_b ] )
    pass
    blat_out_DataFrame_match = blat_out_DataFrame.copy()
    blat_out_DataFrame_match['match_ratio'] = blat_out_DataFrame_match.apply( lambda x: float( '{0:.3f}'.format( (x['query_end'] - x['query_start'])/x['query_size'] ) ), axis = 1 )
    blat_out_DataFrame_match['is_full'] = blat_out_DataFrame_match.apply( lambda x:"Y" if x['match_ratio'] > 0.99 else "N", axis = 1 )
    blat_out_DataFrame_match['is_half'] = blat_out_DataFrame_match.apply( lambda x:"Y" if x['match_ratio'] > 0.50 else "N", axis = 1 )
    blat_out_DataFrame_match_f = blat_out_DataFrame_match.query( 'strand == "+"' )
    blat_out_DataFrame_match_b = blat_out_DataFrame_match.query( 'strand == "-"' )
    if blat_out_DataFrame_match_b.shape[0] > 0:
        blat_out_DataFrame_match_b = blat_out_DataFrame_match_b.assign( new_query_start = lambda x:x['query_size'] - x['query_end'], new_query_end = lambda x:x['query_size'] - x['query_start'] )
        blat_out_DataFrame_match_b = blat_out_DataFrame_match_b.assign(  query_start = lambda x:x['new_query_start'], query_end = lambda x:x['new_query_end'] )
        blat_out_DataFrame_match_b = blat_out_DataFrame_match_b.drop( ["new_query_start", "new_query_end"] , axis=1)
    blat_out_DataFrame_match   = pd.concat( [blat_out_DataFrame_match_f, blat_out_DataFrame_match_b] )
    return( blat_out_DataFrame_match, blat_out_DataFrame_gaps, blat_out_DataFrame_dels )
