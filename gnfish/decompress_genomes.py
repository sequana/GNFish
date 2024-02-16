#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 18:10:29 2021

@author: hectorlorente
"""
import os
import re
import argparse
from gnfish.class_list_files import list_files

import click

import gzip
import shutil




def select_files(directory):
    return list_files(directory)


def decompress_genomes(directory):
    found = False
    for folder in select_files(directory):
        # print(folder[1])
        try:
            for i in range(len(folder[1])):
                file = re.search("(.*).gz", folder[1][i])
                if file:
                    with gzip.open(folder[1][i], "rb") as f_in:
                        with open(file.group(1), "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                            found = True
                            print("Decompressing " + folder[1][i] + " file.")
        except IndexError:
            continue
    if not found:
        print("WARNING. No file decompressed.")


@click.command()
@click.option("--directory", help="sets path to custom folder", type=str)
@click.option("--genomic", help="looks at ./Data/Genomic", is_flag=True)
@click.option("--rna", help="looks at ./Data/Rna", is_flag=True)
@click.option("--protein", help="looks at ./Data/Protein", is_flag=True)
def main(directory, genomic, rna, protein):
    path = os.getcwd()
    data_type = [genomic, rna, protein]
    # Process further based on the options provided

    # def main():
    #     parser = argparse.ArgumentParser(description="Decompress .gz files.")
    #     parser.add_argument("--directory", help="sets path to custom folder", type=str)
    #     parser.add_argument("--genomic", help="looks at ./Data/Genomic", nargs="?", const="genomic", type=str)
    #     parser.add_argument("--rna", help="looks at ./Data/Rna", nargs="?", const="rna", type=str)
    #     parser.add_argument("--protein", help="looks at ./Data/Protein", nargs="?", const="protein", type=str)
    #     args = parser.parse_args()
    #     genomic = args.genomic
    #     rna = args.rna
    #     protein = args.protein
    #     directory = args.directory
    #     path = os.getcwd()
    #     data_type = [genomic, rna, protein]

    if directory is not None:
        decompress_genomes(directory + "/*")
    else:
        for data in data_type:
            if data is not None:
                decompress_genomes(path + "/../Data/" + data.capitalize() + "/*")
