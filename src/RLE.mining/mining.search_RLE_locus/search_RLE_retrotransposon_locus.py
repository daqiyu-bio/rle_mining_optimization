#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------------------#
# Copyright (c) 2026 Beijing Astragenomics Technology Co., Ltd. All rights reserved.        #
# Licensing: see the LICENSE file in the repository root.                                   #
# version : 1.0                                                                             #
#                                                                                           #
# Author:                                                                                   #
#   Daqi Yu <daqi.yu@astragenomics.com>                                                     #
#   Kai Liu <kai.liu@astragenomics.com>                                                     #
#                                                                                           #
# Description:                                                                              #
#   Scan genomes for RLE-associated non-LTR ORF loci, assign RVT families, and export       #
#   grouped locus tables.                                                                   #
# ------------------------------------------------------------------------------------------#

import os
import re
import argparse
import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import traceback

from tools.orfipy_getORF import orfipy_fasta_get_orf_threads
from tools.sequence_get_pfam import sequence_get_pfam
from tools.get_orf_pfam_nonLTR_locus import get_orf_pfam_nonLTR_locus
from tools.get_length_common import get_length_common
from tools.Kmer_analysis import Kmer_generate, Kmer_groupby, Kmer_ratio


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-genome', '--genome_file', required=True, help="CSV describing genome sources (tab-separated columns).")
    parser.add_argument('-hmm', '--RVT_hmm', required=True, help="HMM reference")
    parser.add_argument('-output', '--output_dir', required=True, help="directory to save results")
    parser.add_argument('-log', '--log_dir', required=True, help="directory to save log files")
    args = parser.parse_args()
    return args


def check_dir( temp_dir ):
    os.makedirs(temp_dir,exist_ok=True)
    return temp_dir


def read_fasta(  fasta_path  ):
    with open( fasta_path, 'rt') as fin:
        genome_fna = dict( [(sequence.id, sequence) for sequence in SeqIO.parse(fin, 'fasta')] )
    return( genome_fna )


def get_locus( item_assembly_ID, item_genome_path, item_orf_path, item_orf_pep_path, ref_hmm_path, orf_threads ):
    """
    Predict the ORF and encoded peptide sequences of potential RLE-bearing non-LTR retrotransposons.
    """
    item_fna_orf_pfam       = pd.DataFrame()
    item_fna_orf_pfam_locus = pd.DataFrame()
    item_fna_orf = orfipy_fasta_get_orf_threads( item_genome_path, threads = orf_threads, strand = 'b', minlen=2100, maxlen=6000, starts=['ATG'], stops=['TAA','TAG','TGA'] )
    if item_fna_orf.shape[0] == 0:
        return( item_fna_orf_pfam, item_fna_orf_pfam_locus )
    item_fna_orf['is_X'] = item_fna_orf.apply( lambda x: "Y" if re.search('XXXXXX', x['orf_pep']) else "N", axis = 1)
    item_fna_orf = item_fna_orf.query(  'is_X == "N"' )
    item_fna_orf.insert(0, 'genome_ID', item_assembly_ID )
    item_fna_orf.drop( ['is_X'], axis = 1).to_csv( item_orf_path, sep = "\t" )
    with open( item_orf_pep_path, 'w') as fin:
        for i in range( item_fna_orf.shape[0] ):
            fin.write( "%s\n" % (">" + item_fna_orf['orf_ID'].iloc[i]) )
            fin.write( "%s\n" % (item_fna_orf['orf_pep'].iloc[i]) )
    item_orf_pfam   = sequence_get_pfam( item_orf_pep_path, ref_hmm_path )
    if item_orf_pfam is not None:
        item_fna_orf_pfam = item_fna_orf.merge( item_orf_pfam, left_on = "orf_ID", right_on = "pep_name" ).sort_values( ['chr', 'orf_start'] )
        item_fna_orf_pfam_locus  = get_orf_pfam_nonLTR_locus( item_fna_orf_pfam )
        if item_fna_orf_pfam_locus is not None:
            item_fna_orf_pfam_locus['group'] = Kmer_groupby( item_fna_orf_pfam_locus['orf_pep'], Kmer_generate, Kmer_ratio, Kmer_threshold = 0.6, Kmer_size = 8 )[0]
            item_fna_orf_pfam_locus['group_ID'] = item_fna_orf_pfam_locus.apply( lambda x:x['genome_ID'] + "_nonLTR_" + x['group'], axis=1 )
    return( item_fna_orf_pfam, item_fna_orf_pfam_locus )




def get_genome_locus(item_assembly_ID, item_genome_path, item_orf_path, item_orf_pep_path, ref_hmm_path, item_log_path):
    """
    Get the locus of potential RLE-bearing non-LTR retrotransposons.
    """
    fin = open( item_log_path, 'w' )
    item_log = datetime.datetime.now().strftime( "%F %T" ) + " " + item_assembly_ID + " started."
    fin.write( "%s\n" % (item_log) )
    print( item_log )
    genome_fna_path_size = os.path.getsize( item_genome_path )/( 1024*1024 )
    if  genome_fna_path_size >= 100:
        orf_threads = 8
    else:
        orf_threads = 2
    item_fna_orf_pfam, item_fna_orf_pfam_locus = get_locus( item_assembly_ID, item_genome_path, item_orf_path, item_orf_pep_path, ref_hmm_path, orf_threads )
    item_log = datetime.datetime.now().strftime( "%F %T" ) + " " + item_assembly_ID + " family started."
    fin.write( "%s\n" % (item_log) )
    print( item_log )
    genome_fna_orf_pfam_locus_group = pd.DataFrame()
    if item_fna_orf_pfam_locus is not None and item_fna_orf_pfam_locus.shape[0] > 0:
        item_fna_orf_pfam_locus  = item_fna_orf_pfam_locus.query( 'RLE == RLE' ).copy()
        if item_fna_orf_pfam_locus.shape[0] > 0:
            item_genome_fas   = read_fasta( item_genome_path )
            item_fna_orf_pfam_locus['locus_seq'] = item_fna_orf_pfam_locus.apply( lambda x:str(item_genome_fas[x['chr']].seq)[x['locus_left_end'] - 3000:x['locus_right_end'] + 3000] if x['locus_left_end'] >= 3000 else str(item_genome_fas[x['chr']].seq)[0:x['locus_right_end'] + 3000], axis=1)
            item_fna_orf_pfam_locus['locus_seq'] = item_fna_orf_pfam_locus.apply( lambda x:x['locus_seq'].upper() if x['locus_strand'] == "+" else str( Seq(x['locus_seq']).reverse_complement()).upper(), axis=1)
            item_fna_orf_pfam_locus_groupIDs  = item_fna_orf_pfam_locus.query( 'RLE == RLE' )['group_ID'].drop_duplicates()
            for item in item_fna_orf_pfam_locus_groupIDs:
                item_fna_orf_pfam_locus_group = item_fna_orf_pfam_locus.query( 'group_ID == @item' ).copy()
                item_fna_orf_pfam_locus_group_common = get_length_common( item_fna_orf_pfam_locus_group[ 'orf_pep' ] )[0]
                item_representative = item_fna_orf_pfam_locus_group.iloc[ [ item_fna_orf_pfam_locus_group_common ] ].copy()
                item_representative = item_representative[ ['group_ID', 'genome_ID', 'orf_ID'] ].merge( item_fna_orf_pfam, left_on = ['genome_ID', 'orf_ID'], right_on = ['genome_ID', 'orf_ID'] ).query( 'domain_name == "RVT_1"' )
                item_fna_orf_pfam_locus_group['group_count'] = item_fna_orf_pfam_locus_group.shape[0]
                item_representative = item_representative[['group_ID', 'genome_ID', 'orf_ID']].rename( columns = {"orf_ID":"orf_ID_representative"} )
                item_fna_orf_pfam_locus_group = item_fna_orf_pfam_locus_group.merge( item_representative, left_on = ['genome_ID', 'group_ID'], right_on = ['genome_ID', 'group_ID'] )
                genome_fna_orf_pfam_locus_group = pd.concat( [genome_fna_orf_pfam_locus_group, item_fna_orf_pfam_locus_group] )
    item_log = datetime.datetime.now().strftime( "%F %T" ) + " " + item_assembly_ID + " finished."
    fin.write( "%s\n" % (item_log) )
    print( item_log )
    return( genome_fna_orf_pfam_locus_group )


if __name__ == "__main__":

    args                = get_args()
    genome_path         = args.genome_file
    ref_hmm_path        = args.RVT_hmm
    output_path         = args.output_dir
    log_path            = args.log_dir

    output_path = check_dir(output_path)
    orf_path = check_dir(os.path.join(output_path, 'orf'))
    log_path = check_dir(log_path)

    genomes = pd.read_csv( genome_path, sep="\t" )
    genomes_items = genomes.shape[0]
    outfiles = []
    for i in range( genomes_items ):
        item_assembly_ID  = genomes['Assembly'].iloc[i]
        item_genome_path  = genomes['download_path'].iloc[i]
        item_orf_path     = os.path.join( orf_path, item_assembly_ID + ".orf.csv" )
        item_orf_pep_path = os.path.join( orf_path, item_assembly_ID + ".orf_pep.fasta" )
        item_output_path  = os.path.join( output_path, item_assembly_ID + ".output.csv" )
        item_log_path     = os.path.join( log_path, item_assembly_ID + ".log" )
        
        if os.path.exists( item_log_path ) and os.path.getsize( item_log_path ) > 100 and os.path.exists( item_output_path ):
            outfiles.append(item_output_path)
        else:
            try:
                item_fna_orf_pfam_locus_group = get_genome_locus( item_assembly_ID, item_genome_path, item_orf_path, item_orf_pep_path, ref_hmm_path, item_log_path )
                item_fna_orf_pfam_locus_group.to_csv( item_output_path, sep = "\t" )
                outfiles.append(item_output_path)
            except Exception as error:
                print(traceback.format_exc()) 
                continue
    pd.DataFrame([outfiles]).T.to_csv(f'{output_path}/genome_orf_file.list',index=False,header=None)
