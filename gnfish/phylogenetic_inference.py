#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 16:32:14 2021

@author: hectorlorente
"""

import re
import argparse
import subprocess


def main():
    # Arguments
    parser = argparse.ArgumentParser(
        description="Phylogenetic inference using IQ-Tree."
    )
    parser.add_argument("iqtree_path", help="path to IQ-TREE bin folder", type=str)
    parser.add_argument(
        "input_file",
        help="alignment (FASTA file with prot or nucl aligned seqs)",
        type=str,
    )
    parser.add_argument(
        "--iqtree_parameters",
        help="command to run IQ-Tree, without output declaration, see next parameter. Default TEST, ULTRAFAST BOOSTRAP 1000, ALRT 1000 . See IQ-Tree manual for more options.",
        nargs="?",
        const="-m TEST -B 1000 -alrt 1000",
        type=str,
        default="-m TEST -B 1000 -alrt 1000",
    )
    parser.add_argument(
        "--output_suffix",
        help="suffixes added to input file. Recommended, model and type of support. Default _TEST_UFBS_alrt",
        nargs="?",
        const="_TEST_UFBS_alrt",
        type=str,
        default="_TEST_UFBS_alrt",
    )
    args = parser.parse_args()
    iqtree_path = re.sub("'", "", args.iqtree_path)
    input_file = re.sub("'", "", args.input_file)
    iqtree_parameters = args.iqtree_parameters
    output_suffix = args.output_suffix

    full_iqtree_command = f"{iqtree_path}iqtree2 -s {input_file} {iqtree_parameters} --prefix {input_file}{output_suffix}"
    iqtree_command = full_iqtree_command.split()
    subprocess.run(iqtree_command)
