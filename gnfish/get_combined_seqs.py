#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 13:10:39 2021

@author: hectorlorente
"""


import os
import re
import argparse
from gnfish.utils import list_files
from loguru import logger


def select_files(path):
    return list_files(path)


def get_files(path, pattern):
    files = []
    for folder in select_files(path):
        try:
            found = False
            for i in range(len(folder[1])):
                file = re.search("(.*)" + pattern + "_translated.fas", folder[1][i])
                if file and (not re.search(".gz", file.group())):
                    found = True
                    in_file = re.sub("\\\\", "/", file.group())
                    files.append(in_file)
            if not found:
                for i in range(len(folder[1])):
                    file = re.search("(.*)" + pattern, folder[1][i])
                    if file and (not re.search(".gz", file.group())):
                        found = True
                        in_file = re.sub("\\\\", "/", file.group())
                        files.append(in_file)
            if not found:
                for i in range(len(folder[1])):
                    file = re.search("(.*)RAW_translated.fas", folder[1][i])
                    if file and (not re.search(".gz", file.group())):
                        found = True
                        in_file = re.sub("\\\\", "/", file.group())
                        files.append(in_file)
            if not found:
                for i in range(len(folder[1])):
                    file = re.search("(.*)RAW.fas", folder[1][i])
                    if file and (not re.search(".gz", file.group())):
                        found = True
                        in_file = re.sub("\\\\", "/", file.group())
                        files.append(in_file)
        except IndexError:
            continue
    return files


def read_FASTA_strings(fasta_input_file):
    with open(fasta_input_file) as file:
        return file.read().split(">")[1:]


def read_FASTA_entries(fasta_input_file):
    return [seq.partition("\n") for seq in read_FASTA_strings(fasta_input_file)]


def read_FASTA_sequences(fasta_input_file):
    return [
        [info, seq.replace("\n", "")]
        for info, ignore, seq in read_FASTA_entries(fasta_input_file)
    ]


def combine_all_sequences(
    path, data_type, pattern, output_path, output_name, out_exten
):
    files_lst = get_files(path, pattern)
    with open(output_path + output_name + out_exten, "a") as file:
        for filename in files_lst:
            logger.debug(f"{filename}")
            for line in read_FASTA_sequences(filename):
                file.write(">")
                file.write(line[0])
                file.write("\n")
                file.write(line[1])
                file.write("\n")


def main():
    # Arguments
    parser = argparse.ArgumentParser(
        description="Combines FASTA files in just one matrix"
    )
    parser.add_argument(
        "output_name",
        help="name for the file of the combined matrix",
        nargs="?",
        type=str,
    )
    parser.add_argument(
        "output_path", help="path were to store the combined matrix.", type=str
    )
    parser.add_argument(
        "--genomic",
        help="looks at ./Data/Genomic.",
        nargs="?",
        const="genomic",
        type=str,
    )
    parser.add_argument(
        "--rna", help="looks at ./Data/Rna.", nargs="?", const="rna", type=str
    )
    parser.add_argument(
        "--protein",
        help="looks at ./Data/Protein.",
        nargs="?",
        const="protein",
        type=str,
    )
    parser.add_argument("--directory", help="sets path to custom folder.", type=str)
    parser.add_argument(
        "--pattern",
        help="By default, the pattern for file searching is “final.fas”. As back up to look for “final_translated.fas”, “RAW.fas”, and “RAW_translated.fas”.",
        nargs="?",
        default="final.fas",
        const="final.fas",
        type=str,
    )
    parser.add_argument(
        "--out_exten",
        help='extension of the output file. Default "_all_combined.fas".',
        nargs="?",
        type=str,
        default="_all_combined.fas",
    )
    args = parser.parse_args()
    output_name = args.output_name
    output_path = args.output_path
    genomic = args.genomic
    rna = args.rna
    protein = args.protein
    directory = args.directory
    pattern = args.pattern
    # pattern = re.sub("\.fas", "", pattern)
    out_exten = args.out_exten
    path = os.getcwd()
    data_type_lst = [genomic, rna, protein]
    print(directory)

    if directory is not None:
        for data_type in data_type_lst:
            logger.info(f"Running directory and {data_type} argument.")
            combine_all_sequences(
                directory + "/*",
                data_type,
                pattern,
                output_path,
                output_name,
                out_exten,
            )
    else:
        found = False
        for data_type in data_type_lst:
            if data_type is not None:
                found = True
                logger.info(f"Running {data_type} argument.")
                combine_all_sequences(
                    path + "/../Data/" + data_type.capitalize() + "/*",
                    data_type,
                    pattern,
                    output_path,
                    output_name,
                    out_exten,
                )
        if not found:
            logger.error("You must use genomic, rna or directory arguments.")
