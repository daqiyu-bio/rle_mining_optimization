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
#   Run pyhmmer hmmsearch on a peptide FASTA against a Pfam HMM library.                    #
# ------------------------------------------------------------------------------------------#

import os
import pyhmmer
import pyhmmer.easel
import pyhmmer.plan7
import pandas as pd


def sequence_get_pfam( sequence_fas, hmm_library ):
    aa  = pyhmmer.easel.Alphabet.amino()
    if os.path.exists( sequence_fas ) and os.path.exists( hmm_library ):
        pass
    else:
        print("Input FASTA or HMM file does not exist.")
        return( None )
    with pyhmmer.easel.SequenceFile( sequence_fas, format = 'fasta') as seq_file:
        sequences = list( seq_file )
    sequences = [sequences[i].digitize( aa ) for i in range(len(sequences))]
    pass
    with pyhmmer.plan7.HMMFile( hmm_library ) as hmm_file:
        all_hits = list( pyhmmer.hmmsearch(hmm_file, sequences, cpus=8) )
    pfam_sequence_res = pd.DataFrame()
    for top_hits in all_hits:
        for hit in top_hits:
            pfam_sequence_hit = pd.DataFrame( [hit.name.decode(), hit.best_domain.alignment.target_from, hit.best_domain.alignment.target_to, hit.best_domain.alignment.hmm_name.decode(), hit.best_domain.alignment.hmm_from, hit.best_domain.alignment.hmm_to, hit.best_domain.score, hit.best_domain.pvalue] ).T
            pfam_sequence_res = pd.concat( [pfam_sequence_res, pfam_sequence_hit] )
    pfam_sequence_res_columns = ["pep_name", "target_start", "target_end", "domain_name", "domain_start", "domain_end", "bit_score", "p_value"]
    if pfam_sequence_res.shape[0] > 0:
        pfam_sequence_res.columns = pfam_sequence_res_columns
        return( pfam_sequence_res.sort_values(['pep_name', 'target_start'], ascending = [True, True]) )
    else:
        return( None )

