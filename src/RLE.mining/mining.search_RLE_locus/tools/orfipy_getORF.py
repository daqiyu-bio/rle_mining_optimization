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
#   Predict ORFs via orfipy (per-sequence and whole-genome, with optional multiprocessing). #
# ------------------------------------------------------------------------------------------#

from Bio import SeqIO
import pandas as pd
import multiprocessing as mp
import orfipy_core


def orfipy_sequence_get_orf( sequence, strand = 'b', minlen=600, maxlen=6000, starts=['ATG'], stops=['TAA','TAG','TGA'] ):
    if len( sequence ) <= minlen:
        return( None )
    orfipy_out = orfipy_core.orfs( str( sequence.seq.upper() ), strand=strand, minlen=minlen, maxlen=maxlen, starts=starts, stops=stops )
    orfipy_out = pd.DataFrame( orfipy_out, columns = ['orf_start', 'orf_end', 'orf_strand', 'orf_description'] )
    if orfipy_out.shape[0] == 0:
        return( None )
    orfipy_out.insert(0, "chr", sequence.id )
    orfipy_out['orf_ID']  = range( orfipy_out.shape[0] )
    orfipy_out['orf_ID']  = orfipy_out.apply( lambda x:x['chr'] + "_ORF." + str(x['orf_ID']), axis = 1 )
    orfipy_out['orf_len'] = orfipy_out.apply( lambda x:x['orf_end'] - x['orf_start'], axis = 1)
    orfipy_out['orf_seq'] = orfipy_out.apply( lambda x: str( sequence[ x['orf_start'] : x['orf_end']].seq.upper() ) if x['orf_strand'] == "+" else str( sequence[ x['orf_start'] : x['orf_end']].reverse_complement().seq.upper() ), axis = 1)
    orfipy_out['orf_pep'] = orfipy_out.apply( lambda x: str( sequence[ x['orf_start'] : x['orf_end']].seq.upper().translate() ) if x['orf_strand'] == "+" else str( sequence[ x['orf_start'] : x['orf_end']].reverse_complement().seq.upper().translate() ), axis = 1) 
    orfipy_out['orf_pep_len'] = orfipy_out.apply( lambda x:len( x['orf_pep'] ), axis = 1)
    return( orfipy_out.drop( ['orf_description'], axis = 1 ) )


def orfipy_fasta_get_orf_threads( fasta_path, threads = 8, strand = 'b', minlen=600, maxlen=6000, starts=['ATG'], stops=['TAA','TAG','TGA'] ):
    with open( fasta_path, 'rt') as fin:
        genome_fna = dict( [(sequence.id, sequence) for sequence in SeqIO.parse(fin, 'fasta')] )
    genome_fna_orf = pd.DataFrame() 
    genome_fna_chr = [ i for i in genome_fna.keys() if len(genome_fna[i]) >= minlen ]
    if threads >= mp.cpu_count() / 4:
        threads = int( mp.cpu_count() / 4 )
    if threads >= len( genome_fna_chr ):
        threads = len( genome_fna_chr )
    pool    = mp.Pool( threads )
    results = [ pool.apply_async(orfipy_sequence_get_orf, args=( genome_fna[ item_chr ], strand, minlen, maxlen, starts, stops )) for item_chr in genome_fna_chr ]
    results = [ p.get() for p in results ]
    for item_chr_orf in results:
        genome_fna_orf = pd.concat( [genome_fna_orf, item_chr_orf] )
    return( genome_fna_orf )

