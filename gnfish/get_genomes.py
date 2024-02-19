#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import urllib
from Bio import Entrez
import re
import argparse
import csv

from loguru import logger


# Exceptions
class Error(Exception):
    pass


class ArgumentError(Error):
    def __init__(self, message):
        self.message = message


def read_tsv_file(tsv_file):
    rows = []
    with open(tsv_file) as file:
        tsv_reader = csv.reader(file, delimiter="\t")
        for row in tsv_reader:
            rows.append(row)
        return rows


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
        logger.warning(f"{directory} directory already exits. New data will be added to previous one.")


def get_assembly_data(url, label, species, id_num, path, data_type="genomic", ext=".fna.gz"):
    try:
        link = os.path.join(url, label + "_" + data_type + ext)
        link = re.sub("\\\\", "/", link)
        directory = data_type.capitalize() + "/" + species
        dir_path = re.sub("\\\\", "/", path + "/" + directory)
        create_directory(directory, path)
        urllib.request.urlretrieve(link, "%s/%s" % (dir_path, species + "_" + id_num + "_" + data_type + ext))
        return True
    except:
        dir_path = re.sub("\\\\", "/", path + "/" + directory)
        check_directory(dir_path)
        logger.warning(
            f"No available {data_type} assembly data for {species} {id_num} or unable to download it. Check {species} {id_num} on 'https://www.ncbi.nlm.nih.gov/' assembly database."
        )
        return False


def print_running_info(argument):
    logger.info(f"Running --{argument} argument. Check help or README file for further information.")


def print_already_downloaded_info(data_type, species, id_num):
    logger.info(
        f"{data_type.capitalize()} data for {species} is already downloaded. Check Data/{data_type.capitalize()} directory."
    )


def print_successful_download(data_type, species, id_num):
    logger.info(
        f"{data_type.capitalize()} data for {species}_{id_num} was successfully downloaded. Check Data/{data_type.capitalize()} directory."
    )


def launch_get_assembly_data(url, label, species, id_num, result, exclusive, data_lst, path):
    downloaded = False
    for data_type in data_lst:
        if data_type == "protein":
            if result[2] == "0":
                downloaded = get_assembly_data(url, label, species, id_num, path, data_type, ext=".faa.gz")
            else:
                print_already_downloaded_info(data_type, species, id_num)
            if downloaded:
                result[2] = "1"
                print_successful_download(data_type, species, id_num)
            else:
                if not exclusive:
                    downloaded = get_assembly_data(url, label, species, id_num, path)
                    if downloaded:
                        result[0] = "1"
                        print_successful_download("genomic", species, id_num)
        elif data_type == "genomic":
            if result[0] == "0":
                downloaded = get_assembly_data(url, label, species, id_num, path, data_type)
            else:
                print_already_downloaded_info(data_type, species, id_num)
            if downloaded:
                result[0] = "1"
                print_successful_download(data_type, species, id_num)
        else:
            if result[1] == "0":
                downloaded = get_assembly_data(url, label, species, id_num, path, data_type)
            else:
                print_already_downloaded_info(data_type, species, id_num)
            if downloaded:
                result[1] = "1"
                print_successful_download(data_type, species, id_num)
            else:
                if not exclusive:
                    downloaded = get_assembly_data(url, label, species, id_num, path)
                    if downloaded:
                        result[0] = "1"
                        print_successful_download("genomic", species, id_num)
    return result


def upddate_rows_info(url, label, species, id_num, rows, exclusive, data_lst, path):
    found = False
    for row in rows:
        if row[1] == id_num:
            found = True
            result = launch_get_assembly_data(url, label, species, id_num, row[2:5], exclusive, data_lst, path)
            row[2] = result[0]
            row[3] = result[1]
            row[4] = result[2]
    if not found:
        result = launch_get_assembly_data(url, label, species, id_num, ["0", "0", "0"], exclusive, data_lst, path)
        rows.append([species, id_num, result[0], result[1], result[2]])
    return rows


def get_assembly_summary(id_num, db):
    from Bio import Entrez

    esummary_handle = Entrez.esummary(db=db, id=id_num, report="full")
    esummary_record = Entrez.read(esummary_handle, validate=False)
    return esummary_record


def get_seq_data(term, rows, retmax, db, email, exclusive, data_lst, path):
    Entrez.email = email
    handle = Entrez.esearch(db=db, term=term, retmax=retmax, sort="Significance")
    record = Entrez.read(handle)
    handle.close()
    try:
        error_sentence = record["ErrorList"]["PhraseNotFound"]
        logger.error(
            f"{error_sentence} from {term} not found and can be unpredictable. Check out '{term}' on 'https://www.ncbi.nlm.nih.gov/' assembly database or correct {error_sentence} terms.\n"
        )
    except KeyError:
        ids = record["IdList"]
        for id_num in ids:
            summary = get_assembly_summary(id_num, db)
            species = summary["DocumentSummarySet"]["DocumentSummary"][0]["Organism"]
            species = re.sub(" ", "_", species)
            species = re.search("(.*?_.*?)_.*", species)
            species = species.group(1)
            url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
            if url == "":
                url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]
            if url == "":
                logger.warning(
                    f"No available {db} data for {species}. Check {species} in 'https://www.ncbi.nlm.nih.gov/'"
                )
            else:
                label = os.path.basename(url)
                rows = upddate_rows_info(url, label, species, id_num, rows, exclusive, data_lst, path)
    return rows


def launch_get_genomes(rows, refine, retmax, query_file, db, email, exclusive, data_lst, path):
    new_rows = []
    query_lst = query_file.read().split("\n")
    for query in query_lst:
        if query == "":
            continue
        else:
            term = query + " " + str(refine)
            logger.info(f"Searching for {term} at NCBI Assembly database.")
            new_rows = get_seq_data(term, rows, retmax, db, email, exclusive, data_lst, path)
    return new_rows


def check_file(input_file):
    if os.path.isfile(input_file):
        logger.warning(f"{input_file} file already exists.")
        logger.debug(
            f"New extraction results will be concatenated at the end of {input_file} and the already extant ones will be updated."
        )
        return True
    else:
        return False


def write_extraction_file(tsv_output_file, refine, retmax, query_file, db, email, exclusive, data_lst, path):
    rows = []
    if check_file(tsv_output_file):
        rows = read_tsv_file(tsv_output_file)
    else:
        rows = [["Query_name", "Assembly_ID", "Genomic", "RNA", "Protein"]]
    new_rows = launch_get_genomes(rows, refine, retmax, query_file, db, email, exclusive, data_lst, path)
    # new_rows = new_rows.sort()
    with open(tsv_output_file, "w") as file:
        for new_row in new_rows:
            for col in new_row:
                if col != "":
                    file.write(str(col) + "\t")
            file.write("\n")
    logger.info(
        "Genomes download ended. At Data/downloaded_genomes_log.tsv. You can find information about the downloaded genomes."
    )


def main():
    # Arguments
    parser = argparse.ArgumentParser(
        description="Downloads genome, transcript and protein information from NCBI database."
    )
    parser.add_argument("email", help="mandatory e-mail for NCBI searches", type=str)
    parser.add_argument("query", help="file with the queries", type=argparse.FileType("r"))
    parser.add_argument("--genomic", help="downloads whole genomic data", action="store_true")
    parser.add_argument("--rna", help="downloads protein annotation data", action="store_true")
    parser.add_argument("--protein", help="downloads protein annotation data", action="store_true")
    parser.add_argument(
        "--exclusive", help="download just protein or rna annotation data if available.", action="store_true"
    )
    parser.add_argument(
        "--retmax",
        help="number of NCBI records reported for every query. Default value equal to 200",
        nargs="?",
        const=200,
        type=int,
        default=200,
    )
    parser.add_argument(
        "--refine",
        help='adds filter or field information to all queries. Constant value AND (latest[filter] AND "representative genome"[filter] AND all[filter] NOT anomalous[filter]). Follow constant value structure for your custom refine.',
        nargs="?",
        const='AND (latest[filter] AND "representative genome"[filter] AND all[filter] NOT anomalous[filter])',
        type=str,
        default="",
    )
    args = parser.parse_args()
    query_file = args.query
    genomic = args.genomic
    rna = args.rna
    protein = args.protein
    email = args.email
    exclusive = args.exclusive
    retmax = args.retmax
    refine = args.refine
    data_lst = []
    db = "assembly"
    path = os.getcwd()

    for arg in args.__dict__:
        if args.__dict__[arg]:
            print_running_info(arg)
    create_directory("../Data", path)
    path = path + "/../Data"
    create_directory("Genomic", path)
    create_directory("Rna", path)
    create_directory("Protein", path)

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
    write_extraction_file(
        path + "/downloaded_genomes_log.tsv", refine, retmax, query_file, db, email, exclusive, data_lst, path
    )
