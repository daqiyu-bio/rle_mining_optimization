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
#   Small utilities.                                                                        #
# ------------------------------------------------------------------------------------------#

import os
from subprocess import getstatusoutput
import pandas as pd
from Bio import SeqIO


def run_shell(cmd, quiet=False):
    """
    getstatusoutput:
        - Returns 0 if the command ended successfully;
        - Returns non-zero if the command ended with an error;
        - result contains both stderr and stdout;
    """
    if not quiet:
        print(f'to run: {cmd}')
    Fail, result = getstatusoutput(cmd)
    if not Fail:
        return result
    else:
        print(cmd, 'Execution failed...')
        raise RuntimeError(result)


def mkdir(path):
    """
    Create a directory if it does not exist.
    """
    if os.path.isdir(path):
        print(f'Directory already exists: {path}')
        return path
    else:
        os.makedirs(path)
        print(f'Directory created successfully: {path}')
        return path


def read_fasta(fasta):
    """
    Read a fasta file and return its contents as a dictionary.
    """
    with open(fasta, 'rt') as fin:
        genome_fna = dict([(sequence.id, sequence) for sequence in SeqIO.parse(fin, 'fasta')])
    return genome_fna


def read_psl(psl):
    """
    Read a BLAT PSL result file (*.psl) into a pandas DataFrame.
    """
    blat_out = pd.read_csv(psl, index_col=False, skiprows=5, sep='\s+',
                           names=['match', 'mismatch', 'rep_match', 'N', 'query_gap_count', 'query_gap_base', 'target_gap_count', 'target_gap_base',
                                  'strand', 'query', 'query_size', 'query_start', 'query_end', 'target', 'target_size', 'target_start', 'target_end',
                                  'block_count', 'block_size', "block_query_start", "block_target_start"])
    return blat_out
