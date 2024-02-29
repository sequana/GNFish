#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 13:10:39 2021

@author: hectorlorente
"""


import os
import re
import click
from gnfish.utils import list_files
from loguru import logger


def get_files(path, pattern):
    files = []
    for folder in list_files(path):
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


@click.command()
@click.argument("output_name", type=click.STRING)
@click.argument("output_path", type=click.STRING)
@click.option(
    "--genomic",
    help="looks at ./Data/Genomic.",
    type=click.STRING,
    is_flag=True,
    flag_value="genomic",
)
@click.option(
    "--rna",
    help="looks at ./Data/Rna.",
    type=click.STRING,
    is_flag=True,
    flag_value="rna",
)
@click.option(
    "--protein",
    help="looks at ./Data/Protein.",
    type=click.STRING,
    is_flag=True,
    flag_value="protein",
)
@click.option(
    "--directory",
    help="sets path to custom folder.",
    type=click.Path(exists=True, dir_okay=True),
)
@click.option(
    "--pattern",
    help="By default, the pattern for file searching is “final.fas”. As back up to look for “final_translated.fas”, “RAW.fas”, and “RAW_translated.fas”.",
    type=click.STRING,
    default="final.fas",
    show_default=True,
)
@click.option(
    "--out_exten",
    help="Extension of the output file.",
    type=click.STRING,
    default="_all_combined.fas",
    show_default=True,
)
def main(
    output_name, output_path, genomic, rna, protein, directory, pattern, out_exten
):
    """Combine FASTA files in just one matrix.

    OUTPUT_NAME name for the file of the combined matrix.

    OUTPUT_PATH path were to store the combined matrix.
    """
    path = os.getcwd()
    data_type_lst = [genomic, rna, protein]

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
