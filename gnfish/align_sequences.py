# -*- coding: utf-8 -*-
"""
Program to align sequences using MAFFT tool.
"""

import os
import re
import argparse
import subprocess
from gnfish.utils import list_files

from loguru import logger

import pdb
import pytest

# Arguments
parser = argparse.ArgumentParser(
    description="Aligns files enclosed at Data folder or your custom file. More info at README.md file"
)
parser.add_argument(
    "--genomic",
    help="extracts unique hits from whole genome data stored at Genomic folder or specify genome type if using --directory argument.",
    nargs="?",
    const="genomic",
    type=str,
)
parser.add_argument(
    "--rna",
    help="extracts unique hits from annotated RNA sequences from genomic data stored at Rna folder or specifies genome type if using --directory argument.",
    nargs="?",
    const="rna",
    type=str,
)
parser.add_argument(
    "--protein",
    help="extracts unique hits from annotated protein sequences from genomic data stored at Protein folder or specifies genome type if using --directory argument.",
    nargs="?",
    const="protein",
    type=str,
)
parser.add_argument("--directory", help="sets path to custom folder", type=str)
parser.add_argument(
    "--algorithm",
    help="sets MAFFT algorithm. Default --auto (automatic)",
    nargs="?",
    default="--auto",
    const="--auto",
    type=str,
)
parser.add_argument(
    "--pattern",
    help="sets custom pattern to find files for alignment. Default RAW.fas.",
    nargs="?",
    default="RAW.fas",
    const="RAW.fas",
    type=str,
)
args = parser.parse_args()
genomic = args.genomic
rna = args.rna
protein = args.protein
directory = args.directory
pattern = args.pattern
algorithm = args.algorithm
path = os.getcwd()
data_type_lst = [genomic, rna, protein]


def select_files(path):
    return list_files(path)


def get_files(path):
    files = []
    for folder in select_files(path):
        try:
            found = False
            for i in range(len(folder[1])):
                file = re.search("(.*)" + pattern + ".*", folder[1][i])
                if file and (not re.search(".gz", file.group())):
                    found = True
                    in_file = re.sub("\\\\", "/", file.group())
                    files.append(in_file)
            if not found:
                for i in range(len(folder[1])):
                    file = re.search("(.*)RAW.fas", folder[1][i])
                    if file and (not re.search(".gz", file.group())):
                        in_file = re.sub("\\\\", "/", file.group())
                        files.append(in_file)
        except IndexError:
            continue
    return files


def read_FASTA_strings(fasta_input_file):
    nu_search = re.search("\.fas", fasta_input_file)
    if not nu_search:
        return fasta_input_file.split(">")[1:]
    else:
        with open(fasta_input_file) as file:
            return file.read().split(">")[1:]


def read_FASTA_entries(fasta_input_file):
    return [seq.partition("\n") for seq in read_FASTA_strings(fasta_input_file)]


def read_FASTA_sequences(fasta_input_file):
    return [
        [info, seq.replace("\n", "")]
        for info, ignore, seq in read_FASTA_entries(fasta_input_file)
    ]


if directory is not None:
    found = True
    logger.info("Running directory argument.")
    file_lst = get_files(directory + "/*")
    for infile in file_lst:
        # out = write_mafft_alignment_fasta_file(infile)
        mafft_command = ["/usr/bin/mafft", algorithm, infile]
        output_file = re.search("(.*?)\.fas", infile)
        with open(output_file.group(1) + "_ali.fas", "w") as output:
            subprocess.run(mafft_command, stdout=output, text=True)

        #     file.write(out)
    # else:
    #     found = False
    #     for data_type in data_type_lst:
    #         if data_type is not None:
    #             found = True
    #             logger.info(f"Running data type {data_type} argument.")
    #             file_lst = get_files(path + "/../Data/" + data_type.capitalize() + "/*")
    #             for infile in file_lst:
    #                 out = write_mafft_alignment_fasta_file(infile)
    #                 output_file = re.search("(.*?)\.fas", infile)
    #                 logger.debug(f"{output_file.group(1)}")
    #                 with open(output_file.group(1) + "_ali.fas", "w") as file:
    #                     file.write(out)
    if not found:
        logger.error("You must use genomic, rna or protein arguments.")


def main():
    return ""
