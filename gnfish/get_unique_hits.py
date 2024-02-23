#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 09:56:24 2021

@author: hectorlorente
"""
import os
import sys
import re
import click
import csv
from gnfish.utils import list_files

from loguru import logger


def open_TSV_file(input_file):
    rows = dict()
    with open(input_file) as file:
        tsv_reader = csv.reader(file, delimiter="\t")
        return [row for row in tsv_reader]


def get_files(path, pattern):
    files = []
    found = False
    for folder in list_files(path):
        try:
            for i in range(len(folder[1])):
                file = re.search("(.*)" + pattern, folder[1][i])
                if (
                    file
                    and (not re.search("_unique.tsv", file.group()))
                    and (not re.search("_filtered.tsv", file.group()))
                ):
                    in_file = re.sub("\\\\", "/", file.group())
                    files.append(in_file)
                    found = True

        except IndexError:
            continue
    if not found:
        logger.error(
            f"No file. Try another directory or pattern. Current pattern extension is {pattern}"
        )
        sys.exit()
    else:
        return files


def get_unique_hits(rows):
    new_rows = []
    rows_1 = []
    first = True
    for row in rows:
        if first:
            new_rows.append(row)
            rows_1.append(row[1])
            first = False
        else:
            if not (row[1] in rows_1):
                new_rows.append(row)
                rows_1.append(row[1])
    return new_rows


def generate_output_file(path, pattern):
    tsv_file_lst = get_files(path, pattern)
    for tsv_file in tsv_file_lst:
        new_rows = []
        exten = ""
        logger.debug(f"Extracting hits from {tsv_file}")
        rows = open_TSV_file(tsv_file)
        # if len_treshold:
        #     print('Filtering hits from %s.\n' %tsv_file[0])
        #     exten = '_filtered'
        #     if get_unique:
        #         new_rows = get_unique_hits(rows)
        #         exten = exten + '_unique'
        # elif get_unique:
        new_rows = get_unique_hits(rows)
        exten = "unique.tsv"
        # else:
        #     logger.error('You must use "get_unique" or "filter" arguments.')
        #     sys.exit()
        # print('aqui')
        # print(re.sub(pattern, '_'+ exten, tsv_file))
        with open(re.sub(pattern, "_" + exten, tsv_file), "w") as file:
            for row in new_rows:
                for i in range(len(row)):
                    file.write(row[i] + "\t")
                file.write("\n")


@click.command()
@click.option(
    "--genomic",
    help="Extracts unique hits from genomic data stored at Genomic folder. Use --directory for a custom folder.",
    is_flag=True,
    flag_value="genomic",
    type=str,
)
@click.option(
    "--rna",
    help="Extracts unique hits from rna data stored at Rna folder. Use --directory for a custom folder.",
    is_flag=True,
    flag_value="rna",
    type=str,
)
@click.option(
    "--protein",
    help="Extracts unique hits from protein data stored at Protein folder. Use --directory for a custom folder.",
    is_flag=True,
    flag_value="protein",
    type=str,
)
@click.option("--directory", help="Sets path to a custom folder", type=str)
@click.option(
    "--pattern",
    help='Custom pattern to find Blast output files. Default ".tsv".',
    default=".tsv",
    type=str,
)
def main(genomic, rna, protein, directory, pattern):
    """Gets unique hits from BLAST output files based on genomes IDs."""
    path = os.getcwd()
    data_type_lst = [genomic, rna, protein]

    if directory is not None:
        logger.info("Running directory argument. Blast files must have .tsv extension.")
        generate_output_file(directory + "*", pattern)
    else:
        found = False
        for data_type in data_type_lst:
            if data_type is not None:
                found = True
                logger.info(
                    f"Running {data_type} argument. Blast files must have {pattern} extension."
                )
                generate_output_file(
                    path + "/../Data/" + data_type.capitalize() + "/*", pattern
                )
        if not found:
            logger.error("You must use genomic, rna, protein or directory arguments.")
