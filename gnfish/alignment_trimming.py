#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 21:41:03 2021

@author: hectorlorente
"""

import re
import argparse
import subprocess


def main():
    # Arguments
    parser = argparse.ArgumentParser(description="Alignment trimming using trimAl.")
    parser.add_argument("trimal_path", help="path to trimAl directory", type=str)
    parser.add_argument(
        "input_file",
        help="alignment (FASTA file with prot or nucl aligned seqs)",
        type=str,
    )
    parser.add_argument("output_file", help="output file name.", nargs="?", type=str)
    parser.add_argument(
        "--trimal_parameters",
        help="command to run trimal options. Default -gt 0.1",
        nargs="?",
        const="-gt 0.1",
        type=str,
        default="-gt 0.1",
    )
    args = parser.parse_args()
    trimal_path = re.sub("'", "", args.trimal_path)
    input_file = re.sub("'", "", args.input_file)
    output_file = re.sub("'", "", args.output_file)
    trimal_parameters = args.trimal_parameters

    full_trimal_command = (
        f"{trimal_path}trimal -in {input_file} -out {output_file} {trimal_parameters}"
    )
    trimal_command = full_trimal_command.split()
    subprocess.run(trimal_command)
