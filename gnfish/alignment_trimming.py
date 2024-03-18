#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 21:41:03 2021

@author: hectorlorente
"""

import click
import subprocess


@click.command()
@click.argument("trimal_path", type=click.STRING)
@click.argument(
    "input_file",
    type=click.STRING,
)
@click.argument("output_file", type=click.STRING)
@click.option(
    "--trimal_parameters",
    help="Command to run trimal options.",
    type=click.STRING,
    default="-gt 0.1",
    show_default=True,
)
def main(trimal_path, input_file, output_file, trimal_parameters):
    """Trim an alignment using trimAl.

    TRIMAL_PATH Path to trimAl directory.

    INPUT_FILE Alignment (FASTA file with prot or nucl aligned seqs.

    OUTPUT_FILE Output file name.

    """
    full_trimal_command = (
        f"{trimal_path}trimal -in {input_file} -out {output_file} {trimal_parameters}"
    )
    trimal_command = full_trimal_command.split()
    subprocess.run(trimal_command)
