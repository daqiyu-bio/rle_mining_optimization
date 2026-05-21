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
#   K-mer similarity grouping for peptide clustering.                                       #
# ------------------------------------------------------------------------------------------#

def Kmer_generate( string, Kmer_size = 8):
    if len(string) >= Kmer_size:
        Kmer_set = set( string[i:(i + Kmer_size)] for i in range( len(string) - Kmer_size )  )
    else:
        Kmer_set = set( string )
    return( Kmer_set )


def Kmer_ratio( Kmer_set1, Kmer_set2 ):
    return( max( len( Kmer_set1.intersection( Kmer_set2) )/len( Kmer_set1 ), len( Kmer_set1.intersection( Kmer_set2 ) )/len( Kmer_set2 ) ) )


def Kmer_groupby( string_list, Kmer_generate, Kmer_ratio, Kmer_threshold = 0.5, Kmer_size = 8 ):
    Kmer_list = [ Kmer_generate(i, Kmer_size) for i in string_list]
    Kmer_list_len = len( Kmer_list )
    Kmer_groups = []
    Kmer_unions = []
    Kmer_ratios = []
    Group_lists  = []
    for i in range( Kmer_list_len ):
        if i == 0:
            Group_lists.append( i )
            Kmer_groups.append( [i] )
            Kmer_unions.append( Kmer_list[i] )
            Kmer_ratios.append( 0 )
        else:
            Group_mark = 0
            for j in range( len(Kmer_unions) ):
                s = Kmer_ratio( Kmer_list[i], Kmer_unions[j] )
                if s > Kmer_threshold:
                    Group_lists.append( j )
                    Kmer_groups[j].append( i )
                    Kmer_unions[j].update( Kmer_list[i] )
                    Group_mark = 1
                    Kmer_ratios.append( s )
                    break
            if Group_mark == 0:
                Group_lists.append( len(Kmer_unions) )
                Kmer_groups.append( [i] )
                Kmer_unions.append( Kmer_list[i] )
                Kmer_ratios.append( 0 )
    Group_lists = [ "g" + str(i) for i in Group_lists ]
    return( Group_lists, Kmer_ratios )
