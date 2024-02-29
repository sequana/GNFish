#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import urllib
from Bio import Entrez
import re
import click
import pandas as pd
from loguru import logger

import pdb


# Exceptions
class Error(Exception):
    pass


class ArgumentError(Error):
    def __init__(self, message):
        self.message = message


def check_directory(directory):
    directory = directory
    directory_content = os.listdir(directory)
    if len(directory_content) == 0:
        os.rmdir(directory)


def create_directory(directory, path):
    directory = path + "/" + directory
    if not os.path.isdir(directory):
        os.mkdir(directory)
    else:
        logger.warning(
            f"{directory} directory already exits. New data will be added to previous one."
        )


def get_assembly_data(url, label, species, id_num, path, data_type, ext=".fna.gz"):
    try:
        if data_type == "protein":
            ext = ".faa.gz"
        link = os.path.join(url, label + "_" + data_type + ext)
        link = re.sub("\\\\", "/", link)
        directory = data_type.capitalize() + "/" + species
        dir_path = re.sub("\\\\", "/", path + "/" + directory)
        create_directory(directory, path)
        urllib.request.urlretrieve(
            link, "%s/%s" % (dir_path, species + "_" + id_num + "_" + data_type + ext)
        )
        return True
    except:
        dir_path = re.sub("\\\\", "/", path + "/" + directory)
        check_directory(dir_path)
        logger.warning(
            f"No available {data_type} assembly data for {species} {id_num} or unable to download it. Check {species} {id_num} on 'https://www.ncbi.nlm.nih.gov/' assembly database."
        )
        return False


def print_running_info(argument):
    logger.info(
        f"Running --{argument} argument. Check help or README file for further information."
    )


def print_already_downloaded_info(data_type, species, id_num):
    logger.info(
        f"{data_type.capitalize()} data for {species} {id_num} is already downloaded. Check Data/{data_type.capitalize()} directory."
    )


def print_successful_download(data_type, species, id_num):
    logger.info(
        f"{data_type.capitalize()} data for {species}_{id_num} was successfully downloaded. Check Data/{data_type.capitalize()} directory."
    )


def launch_get_assembly_data(
    url, label, species, id_num, log_file_data, exclusive, data_lst, path
):
    downloaded = False
    for data_type in data_lst:
        resultado = log_file_data.loc[
            log_file_data["Assembly_ID"] == id_num, data_type.capitalize()
        ]
        if int(resultado.iloc[0]) == 0:
            downloaded = get_assembly_data(url, label, species, id_num, path, data_type)
        else:
            print_already_downloaded_info(data_type, species, id_num)
        if (
            not exclusive
            and data_type != "genomic"
            and not downloaded
            and log_file_data.loc[log_file_data["Assembly_ID"] == id_num, "Genomic"]
            == "0"
        ):
            data_type = "genomic"
            downloaded = get_assembly_data(url, label, species, id_num, path, data_type)
        if downloaded:
            log_file_data.loc[
                log_file_data["Assembly_ID"] == id_num, data_type.capitalize()
            ] = "1"
            print_successful_download("genomic", species, id_num)
    return log_file_data


def get_assembly_summary(id_num, db):
    esummary_handle = Entrez.esummary(db=db, id=id_num, report="full")
    esummary_record = Entrez.read(esummary_handle, validate=False)
    return esummary_record


def check_file(input_file):
    if os.path.isfile(input_file):
        logger.warning(f"{input_file} file already exists.")
        logger.debug(
            f"New extraction results will be concatenated at the end of {input_file} and the already extant ones will be updated."
        )
        return True
    else:
        return False


def extract_url(summary):
    url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
    if not url:
        url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]
    return url


def get_record_entrez(db, email, final_query, retmax):
    Entrez.email = email
    handle = Entrez.esearch(db=db, term=final_query, retmax=retmax, sort="Significance")
    record = Entrez.read(handle)
    handle.close()
    return record


def get_species(summary):
    species = summary["DocumentSummarySet"]["DocumentSummary"][0]["Organism"]
    species = re.sub(" ", "_", species)
    species = re.search("(.*?_.*?)_.*", species)
    species = species.group(1)
    return species


def manage_create_directory(path):
    create_directory("/Data", path)
    path = path + "/Data"
    create_directory("Genomic", path)
    create_directory("Rna", path)
    create_directory("Protein", path)
    return path


@click.command()
@click.argument("working_dir", type=click.STRING, required=True)
@click.argument(
    "email",
    type=click.STRING,
    required=True,
)
@click.argument(
    "query",
    type=click.File("r"),
    required=True,
)
@click.option("--genomic", is_flag=True, help="Downloads whole genomic data.")
@click.option("--rna", is_flag=True, help="Downloads rna annotation data.")
@click.option("--protein", is_flag=True, help="Downloads protein annotation data.")
@click.option(
    "--exclusive",
    is_flag=True,
    help="Download just protein or RNA annotation data if available.",
)
@click.option(
    "--retmax",
    type=click.INT,
    default=200,
    help="Number of NCBI records reported for every query.",
    show_default=True,
)
@click.option(
    "--refine",
    type=click.STRING,
    default="",
    show_default=True,
    is_flag=False,
    flag_value="AND (latest[filter] AND 'representative genome'[filter] AND all[filter] NOT anomalous[filter])",
    help='Adds filter or field information to all queries. Flagq value "AND (latest[filter] AND "representative genome"[filter] AND all[filter] NOT anomalous[filter])". Follow constant value structure for your custom refine.',
)
def main(email, query, genomic, rna, protein, exclusive, retmax, refine):
    """Download genome, transcript, and protein data from NCBI database.

        EMAIL your email for connecting to NCBI database.

    QUERY path to the file with your queries.
    """
    params_to_log = {key: value for key, value in locals().items() if value is not None}
    print_running_info(params_to_log)
    data_lst = []
    db = "assembly"
    path = manage_create_directory(os.getcwd())
    tsv_output_file = path + "/downloaded_genomes_log.tsv"
    rows = []
    query_lst = query.read().split("\n")

    # Checks data type
    if genomic == rna == protein is False:
        data_lst.append("genomic")
    else:
        if rna:
            data_lst.append("rna")
        if protein:
            data_lst.append("protein")
        if genomic:
            data_lst.append("genomic")

    # Runs the program
    if check_file(log_file):
        log_file_data = pd.read_csv(log_file, sep="\t")
        log_file_data["Assembly_ID"] = log_file_data["Assembly_ID"].astype(str)
    else:
        log_file_data = pd.DataFrame(
            columns=["Species", "Assembly_ID", "Genomic", "RNA", "Protein"]
        )
    for query in query_lst:
        if query == "":
            continue
        final_query = query + " " + str(refine)
        logger.info(f"Searching for {final_query} at NCBI Assembly database.")
        record = get_record_entrez(db, email, final_query, retmax)
        try:
            error_sentence = record["ErrorList"]["PhraseNotFound"]
            logger.error(
                f"{error_sentence} from {final_query} not found and can be unpredictable. Check out '{final_query}' on 'https://www.ncbi.nlm.nih.gov/' assembly database or correct {error_sentence} terms.\n"
            )
        except KeyError:
            ids = record["IdList"]
            for id_num in ids:
                summary = get_assembly_summary(id_num, db)
                species = get_species(summary)
                url = extract_url(summary)
                if not url:
                    logger.warning(
                        f"No available {db} data for {species}. Check {species} in 'https://www.ncbi.nlm.nih.gov/'"
                    )
                else:
                    if id_num not in log_file_data["Assembly_ID"].values:
                        new_id_num = {
                            "Species": species,
                            "Assembly_ID": id_num,
                            "Genomic": 0,
                            "RNA": 0,
                            "Protein": 0,
                        }
                        log_file_data.loc[len(log_file_data)] = new_id_num
                    label = os.path.basename(url)
                    log_file_data = launch_get_assembly_data(
                        url,
                        label,
                        species,
                        id_num,
                        log_file_data,
                        exclusive,
                        data_lst,
                        data_dir,
                    )
    log_file_data.to_csv(log_file, sep="\t", index=False)
    logger.info(
        "Genomes download ended. At Data/downloaded_genomes_log.tsv. You can find information about the downloaded genomes."
    )
