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
import subprocess
from gnfish.utils import list_files
from loguru import logger


def get_files(path):
    files = set()
    out_files = {}
    found = False
    for folder in list_files(path):
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
        logger.error(
            'No decompressed genome file found. Decompress your ".gz" files using "decompress_genomes.py".'
        )
        sys.exit()
    else:
        return files


def makeblastdb(db_file, db_type, blast_path):
    if db_type == "protein":
        blastCommand = [blast_path + "makeblastdb", "-dbtype", "prot", "-in", db_file]
    else:
        blastCommand = [blast_path + "makeblastdb", "-dbtype", "nucl", "-in", db_file]
    subprocess.run(blastCommand)


def blast(path, db_type, blast_path, query_type, query_file, out_exten, outfmt, evalue):
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
        makeblastdb(db_file, db_type, blast_path)
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


@click.command()
@click.argument("blast_path", type=str, required=True)
@click.argument("query_file", type=str, required=True)
@click.argument("query_type", type=str, required=True)
@click.option(
    "--genomic",
    type=str,
    help="blast against genomic data stored at Genomic folder. Use directory for custom folder.",
    is_flag=False,
    flag_value="genomic",
)
@click.option(
    "--rna",
    type=str,
    help="blast against rna data stored at Rna folder. Use directory for custom folder.",
    is_flag=False,
    flag_value="rna",
)
@click.option(
    "--protein",
    type=str,
    help="blast against protein data stored at Protein folder. Use directory for custom folder.",
    is_flag=False,
    flag_value="protein",
)
@click.option(
    "--directory",
    type=str,
    help="sets path to custom folder (Genome files must have .fna (genomic and rna data) or .faa (protein data) extensions)",
)
@click.option(
    "--evalue",
    type=str,
    help="E-value threshold for search. Default value equal to 1e-10",
    default="1e-10",
)
@click.option(
    "--outfmt", type=str, help="Output alignment options. Default value 6", default="6"
)
@click.option(
    "--out_exten",
    type=str,
    help='Extension of the BLAST output file. Default "_out.tsv".',
    default="_out.tsv",
)
def main(
    blast_path,
    query_file,
    query_type,
    genomic,
    rna,
    protein,
    directory,
    evalue,
    outfmt,
    out_exten,
):
    """BLAST against genome .fna and .faa files.

    BLAST_PATH path to ncbi-blast directory.

    QUERY_FILE path to your query file with your reference sequences (Fasta file with prot or nucl seqs).

    QUERY_TYPE data type of query sequences. prot or nucl.
    """
    path = os.getcwd()
    db_type_lst = [genomic, rna, protein]

    if re.search("prot|nucl", query_type):
        if directory is not None:
            found = False
            for db_type in db_type_lst:
                if db_type is not None:
                    found = True
                    logger.info(
                        f"Running directory and {db_type} arguments. Genome files must have .fna (genomic and rna data) or .faa (protein data) extensions."
                    )
                    blast(
                        directory + "/*",
                        db_type,
                        blast_path,
                        query_type,
                        query_file,
                        out_exten,
                        outfmt,
                        evalue,
                    )
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
                    blast(
                        path + "/../Data/" + db_type.capitalize() + "/*",
                        db_type,
                        blast_path,
                        query_type,
                        query_file,
                        out_exten,
                        outfmt,
                        evalue,
                    )
            if not found:
                logger.error("You must use genomic, rna or protein arguments.")
    else:
        logger.error('db_type argument must be "prot" or "nucl".')
