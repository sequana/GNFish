#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 16:32:14 2021

@author: hectorlorente
"""

import click
import subprocess


@click.command()
@click.argument("iqtree_path", type=click.STRING)
@click.argument("input_file", type=click.STRING)
@click.option(
    "--iqtree_parameters",
    type=click.STRING,
    help="command to run IQ-Tree, without output declaration, see next parameter. Default TEST, ULTRAFAST BOOSTRAP 1000, ALRT 1000 . See IQ-Tree manual for more options.",
    default="-m TEST -B 1000 -alrt 1000",
)
@click.option(
    "--output_suffix",
    type=click.STRING,
    help="suffixes added to input file. Recommended, model and type of support. Default _TEST_UFBS_alrt",
    default="_TEST_UFBS_alrt",
)
def main(iqtree_path, input_file, iqtree_parameters, output_suffix):
    """Conduct phylogenetic inference using IQ-Tree.

        IQTREE_PATH path to IQ-TREE bin folder.

    INPUT_FILE alignment (FASTA file with prot or nucl aligned seqs.
    """
    full_iqtree_command = f"{iqtree_path}iqtree2 -s {input_file} {iqtree_parameters} --prefix {input_file}{output_suffix}"
    iqtree_command = full_iqtree_command.split()
    subprocess.run(iqtree_command)
