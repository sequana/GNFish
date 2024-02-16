#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 09:56:24 2021

@author: hectorlorente
"""
import os
import sys
import urllib
import re
import argparse
import csv
import subprocess
from class_list_files import list_files
from loguru import logger


##Arguments
parser = argparse.ArgumentParser(description="BLAST against genome .fna and .faa files.")
parser.add_argument("blast_path", help="path to ncbi-blast directory", type=str)
parser.add_argument(
    "query_file",
    help="path to your query file with your reference sequences (Fasta file with prot or nucl seqs)",
    type=str,
)
parser.add_argument("query_type", help="data type of query sequences. prot or nucl", type=str)
parser.add_argument(
    "--genomic",
    help="blast against genomic data stored at Genomic folder. Use directory for custom folder.",
    nargs="?",
    const="genomic",
    type=str,
)
parser.add_argument(
    "--rna",
    help="blast against rna data stored at Rna folder. Use directory for custom folder.",
    nargs="?",
    const="rna",
    type=str,
)
parser.add_argument(
    "--protein",
    help="blast against protein data stored at Protein folder. Use directory for custom folder.",
    nargs="?",
    const="protein",
    type=str,
)
parser.add_argument(
    "--directory",
    help="sets path to custom folder (Genome files must have .fna (genomic and rna data) or .faa (protein data) extensions)",
    type=str,
)
parser.add_argument(
    "--evalue", help="E-value threshold for search. Default value equal to 1e-10", nargs="?", type=str, default="1e-10"
)
parser.add_argument("--outfmt", help="Output alignment options. Default value 6", nargs="?", type=str, default="6")
parser.add_argument(
    "--out_exten",
    help='Extension of the BLAST output file. Default "_out.tsv".',
    nargs="?",
    type=str,
    default="_out.tsv",
)
args = parser.parse_args()
blast_path = re.sub("'", "", args.blast_path)
query_file = re.sub("'", "", args.query_file)
query_type = args.query_type
genomic = args.genomic
rna = args.rna
protein = args.protein
directory = args.directory
evalue = args.evalue
outfmt = args.outfmt
out_exten = args.out_exten
path = os.getcwd()
db_type_lst = [genomic, rna, protein]


def select_files(path):
    return list_files.list_files_method(list_files, path)


def get_files(path):
    files = set()
    out_files = {}
    found = False
    for folder in select_files(path):
        try:
            for i in range(len(folder[1])):
                file = re.search("(((.*)\.f[a,n]a).*)", folder[1][i])
                if file:
                    in_file = re.sub("\\\\", "/", file.group(2))
                    files.add(in_file)
                    found = True
        except IndexError:
            continue
    if not found:
        logger.error('No decompressed genome file found. Decompress your ".gz" files using "decompress_genomes.py".')
        sys.exit()
    else:
        return files


def makeblastdb(db_file, db_type):
    if db_type == "protein":
        blastCommand = [blast_path + "makeblastdb", "-dbtype", "prot", "-in", db_file]
    else:
        blastCommand = [blast_path + "makeblastdb", "-dbtype", "nucl", "-in", db_file]
    subprocess.run(blastCommand)


def blast(path, db_type):
    db_files = get_files(path)
    if db_type == "protein":
        if re.search("prot", query_type):
            blastProgram = "blastp"
        else:
            blastProgram = "blastx"
    else:
        if re.search("prot", query_type):
            blastProgram = "tblastn"
        else:
            blastProgram = "blastn"
    for db_file in db_files:
        makeblastdb(db_file, db_type)
        blastCommand = [
            blast_path + blastProgram,
            "-query",
            query_file,
            "-db",
            db_file,
            "-out",
            re.sub("\.f.*", out_exten, db_file),
            "-outfmt",
            outfmt,
            "-evalue",
            evalue,
        ]
        subprocess.run(blastCommand)


if re.search("prot|nucl", query_type):
    if directory is not None:
        found = False
        for db_type in db_type_lst:
            if db_type is not None:
                found = True
                logger.info(
                    f"Running directory and {db_type} arguments. Genome files must have .fna (genomic and rna data) or .faa (protein data) extensions."
                )
                blast(directory + "/*", db_type)
        if not found:
            logger.error("You must use genomic, rna or protein arguments.")
    else:
        found = False
        for db_type in db_type_lst:
            if db_type is not None:
                found = True
                logger.info(
                    "Running {db_type} argument. Genome files must have .fna (genomic and rna data) or .faa (protein data) extensions."
                )
                blast(path + "/../Data/" + db_type.capitalize() + "/*", db_type)
        if not found:
            logger.error("You must use genomic, rna or protein arguments.")
else:
    logger.error('db_type argument must be "prot" or "nucl".')
