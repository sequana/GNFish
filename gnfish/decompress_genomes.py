#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 18:10:29 2021

@author: hectorlorente
"""
import os
import re
from gnfish.utils import list_files

import click

import gzip
import shutil

from loguru import logger


def decompress_genomes(directory):
    found = False
    for folder in list_files(directory):
        try:
            for i in range(len(folder[1])):
                file = re.search("(.*).gz", folder[1][i])
                if file:
                    with gzip.open(folder[1][i], "rb") as f_in:
                        with open(file.group(1), "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                            found = True
                            logger.info(f"Decompressing {folder[1][i]} file.")
        except IndexError:
            continue
    if not found:
        logger.warning("No file decompressed.")


@click.command()
@click.option("--directory", help="sets path to custom folder", type=str)
@click.option("--genomic", help="looks at ./Data/Genomic", is_flag=True)
@click.option("--rna", help="looks at ./Data/Rna", is_flag=True)
@click.option("--protein", help="looks at ./Data/Protein", is_flag=True)
def main(directory, genomic, rna, protein):
    path = os.getcwd()
    data_type = [genomic, rna, protein]

    if directory is not None:
        decompress_genomes(directory + "/*")
    else:
        for data in data_type:
            if data is not None:
                decompress_genomes(path + "/../Data/" + data.capitalize() + "/*")
