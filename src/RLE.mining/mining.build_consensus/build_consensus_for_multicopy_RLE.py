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
#   Parallel consensus builder for multicopy RLE locus MSAs (reads a list of group_ID       #
#   and MSA path pairs).                                                                    #
# ------------------------------------------------------------------------------------------#

import argparse
import os
import traceback
from multiprocessing import Pool

import pandas as pd
from Bio import Align
from Bio import SeqIO

from align_consensus_parse import align_consensus_parse_RNA


def read_fasta(fasta):
    with open(fasta, 'rt') as fin:
        genome_fna = dict([(sequence.id, sequence) for sequence in SeqIO.parse(fin, 'fasta')])
    return genome_fna


def get_args():
    parser = argparse.ArgumentParser(
        description="Build consensus from per-group locus MSA files."
    )
    parser.add_argument(
        "-msa-list",
        "--msa-list",
        dest="msa_list",
        required=True,
        help="Two columns per line: group_ID then MSA path (tab or spaces; use tab if path contains spaces). Lines starting with # and empty lines are skipped.",
    )
    parser.add_argument(
        "-output",
        "--output",
        dest="output",
        required=True,
        help="Directory to write <group_ID>.consensus.csv / .warning / .error",
    )
    parser.add_argument(
        "-conserve",
        "--conserve",
        dest="conserve",
        action="store_true",
        help="Pass conserve=True into align_consensus_parse_RNA (default: False).",
    )
    parser.add_argument(
        "-mismatch-size",
        "--mismatch-size",
        dest="mismatch_size",
        type=int,
        default=40,
        help="mismatch_size for align_consensus_parse_RNA (default: 40).",
    )
    parser.add_argument(
        "-threads",
        "--threads",
        dest="threads",
        type=int,
        default=60,
        help="Worker processes for Pool (default: 60).",
    )
    return parser.parse_args()


def read_msa_list(path):
    """Each line: group_ID <sep> MSA_file. First column is group_ID; remainder of line is path."""
    rows = []
    with open(path, "rt") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) < 2:
                raise ValueError(
                    f"Expected two columns (group_ID and MSA path), got {len(parts)} field(s): {line!r}"
                )
            rows.append((parts[0], parts[1]))
    return rows


def make_aligner_strict():
    aligner_strict = Align.PairwiseAligner()
    aligner_strict.mode = "local"
    aligner_strict.match_score = 2
    aligner_strict.mismatch_score = -12
    aligner_strict.open_gap_score = -12
    aligner_strict.extend_gap_score = -1
    return aligner_strict


def trim(group_ID, MSA_Path, aligner_strict, conserve, mismatch_size, outpath):
    MSA_raw_size = len(read_fasta(MSA_Path))
    try:
        ret = align_consensus_parse_RNA(
            MSA_Path,
            aligner_strict,
            conserve=conserve,
            mismatch_size=mismatch_size,
        )
        if ret is not None:
            align_consensus_choice, align_consensus_choice_seq, align_consensus_msa_start, align_consensus_msa_end = list(
                ret
            )
            Left_choice = list(align_consensus_choice[0])
            Right_choice = list(align_consensus_choice[1])
            df = pd.DataFrame(
                [
                    [
                        group_ID,
                        MSA_Path,
                        MSA_raw_size,
                        Left_choice,
                        Right_choice,
                        align_consensus_msa_start,
                        align_consensus_msa_end,
                        align_consensus_choice_seq,
                    ]
                ],
                columns=[
                    "group_ID",
                    "MSA_Path",
                    "MSA_raw_size",
                    "Left_choice",
                    "Right_choice",
                    "align_consensus_msa_start",
                    "align_consensus_msa_end",
                    "align_consensus_choice_seq",
                ],
            )

            df.to_csv(f"{outpath}/{group_ID}.consensus.csv", sep="\t", index=False)
        else:
            print(f"{group_ID} warning")
            os.system(
                f"echo can not build consensus. >  {outpath}/{group_ID}.consensus.warning"
            )
            return 0
    except Exception:
        print(f"{group_ID} error")
        error_msg = traceback.format_exc()
        os.system(
            f'echo can not build consensus. "{error_msg}" >  {outpath}/{group_ID}.consensus.error'
        )
        return 0


def main():
    args = get_args()
    msa_rows = read_msa_list(args.msa_list)
    outpath = args.output
    aligner_strict = make_aligner_strict()

    pool = Pool(args.threads)
    for group_ID, MSA_Path in msa_rows:
        pool.apply_async(
            trim,
            args=(
                group_ID,
                MSA_Path,
                aligner_strict,
                args.conserve,
                args.mismatch_size,
                outpath,
            ),
        )
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
