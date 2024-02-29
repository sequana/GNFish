import os
import re
import csv
import click
from gnfish.utils import list_files
from Bio import Entrez
from Bio import SeqIO
from loguru import logger


def get_files(path, pattern):
    files = []
    for folder in list_files(path):
        try:
            for i in range(len(folder[1])):
                file = re.search("(.*)" + pattern + ".*", folder[1][i])
                if (
                    file
                    and (not re.search("\.gz", file.group()))
                    and (not re.search("\.f[a,n]a\.[n,p].*$", file.group()))
                ):
                    in_file = re.sub("\\\\", "/", file.group())
                    files.append(in_file)
        except IndexError:
            continue
    return files


def open_CSV_file(csv_input_file, blast_pattern):
    try:
        with open(csv_input_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter="\t")
            return [[row[0], row[1], row[6], row[7]] for row in csv_reader]
    except FileNotFoundError:
        logger.warning(
            f"{csv_input_file} file not found. Could be due to a problem with '{str(blast_pattern)}' pattern. Check the file or change --blast_pattern parameter."
        )


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


def read_FASTA_strings_2(fasta_input_file):
    return fasta_input_file.split(">")[1:]


def read_FASTA_entries_2(fasta_input_file):
    return [seq.partition("\n") for seq in read_FASTA_strings_2(fasta_input_file)]


def read_FASTA_sequences_2(fasta_input_file):
    return [
        [info, seq.replace("\n", "")]
        for info, ignore, seq in read_FASTA_entries_2(fasta_input_file)
    ]


def reverse_complement(DNA_string):
    DNA_string = DNA_string.upper()
    dic_complement = {"A": "T", "G": "C", "C": "G", "T": "A"}
    reverse_string = str()
    for str_base in DNA_string:
        try:
            reverse_string += dic_complement[str_base]
        except:
            reverse_string += str_base
    return reverse_string[::-1]


def get_seqs_data(data_type, id_num, email):
    Entrez.email = email
    if not re.search("protein", data_type):
        handle = Entrez.efetch(db="protein", id=id_num, rettype="gp", retmode="text")
        record = SeqIO.read(handle, "genbank")
        id_num = record.annotations.get("db_source")
        handle = Entrez.efetch(
            db="nucleotide", id=id_num, rettype="fasta", retmode="text"
        )
        record = handle.read()
    else:
        handle = Entrez.efetch(db="protein", id=id_num, rettype="fasta", retmode="text")
        record = handle.read()
    return record


def crossvalidate_sequences(
    blast_uni_lines_lst,
    blast_lines_lst,
    genome_file,
    data_type,
    blast_uni_file,
    blast_file,
    query_seqs,
    query_seqs_num,
    in_len,
    email,
):
    lines = []
    fasta_lines_lst = read_FASTA_sequences(genome_file)
    length = len(fasta_lines_lst)
    for row in blast_uni_lines_lst:
        ali_lines = []
        row_found = False
        for line in fasta_lines_lst:
            row_search = re.search(row[1], line[0])
            mod_row = re.search("(.*?[0-9])_[0-9]", row[1])
            if row_search or (mod_row and re.search(mod_row.group(1), line[0])):
                row_found = True
                line_0 = re.sub(".*\.[0-9]", row[1], line[0])
                line_1 = line[1]
                if re.search("genomic", data_type):
                    if row[2] < row[3]:
                        lower_bound = int(row[2]) - in_len
                        upper_bound = int(row[3]) + in_len
                        if lower_bound < 0:
                            line_1 = line[1][0:upper_bound]
                        else:
                            line_1 = line[1][lower_bound:upper_bound]
                    else:
                        lower_bound = int(row[3]) - in_len
                        upper_bound = int(row[2]) + in_len
                        if lower_bound < 0:
                            line_1 = reverse_complement(line[1][0:upper_bound])
                        else:
                            line_1 = reverse_complement(
                                line[1][lower_bound:upper_bound]
                            )
                ali_lines.append([line_0, line_1])
                if query_seqs:
                    count = 0
                    bl_row_found = False
                    blast_lines_lst
                    for bl_row in blast_lines_lst:
                        if re.search(row[1], bl_row[1]) and count < query_seqs_num:
                            count += 1
                            bl_row_found = True
                            db_seq = get_seqs_data(data_type, bl_row[0], email)
                            db_seq = read_FASTA_sequences_2(db_seq)
                            ali_lines.append(db_seq)
                    if not bl_row_found:
                        logger.warning(f"Unable to find '{row[1]}' in '{blast_file}'.")
                break
        if not row_found:
            logger.warning('. Unable to find "%s" in "%s".' % (row[1], genome_file))
        lines.append((ali_lines, row[1]))
    return lines


def generate_output_FASTA_file(
    path,
    data_type,
    blast_pattern,
    genome_pattern,
    query_seqs,
    query_seqs_num,
    in_len,
    email,
):
    count = 0
    blast_uni_files_lst = get_files(path, blast_pattern)
    genome_files_lst = get_files(path, genome_pattern)
    for blast_uni_file in blast_uni_files_lst:
        found = True
        ID = re.search("[A-Z]+[a-z]*_[a-z]*?[a-z]*?_([0-9].*?)_.*", blast_uni_file)
        if not ID:
            ID = re.search("(GC[A,F]_[0-9].*?)_", blast_uni_file)
        if ID:
            blast_uni_lines_lst = open_CSV_file(blast_uni_file, blast_pattern)
            if blast_uni_lines_lst is not None:
                found = False
                blast_file = re.search("(.*?)_unique.tsv", blast_uni_file)
                blast_file = blast_file.group(1) + ".tsv"
                blast_lines_lst = open_CSV_file(blast_file, blast_pattern)
                for genome_file in genome_files_lst:
                    if re.search(ID.group(1), genome_file):
                        found = True
                        output_lines_lst = crossvalidate_sequences(
                            blast_uni_lines_lst,
                            blast_lines_lst,
                            genome_file,
                            data_type,
                            blast_uni_file,
                            blast_file,
                            query_seqs,
                            query_seqs_num,
                            in_len,
                            email,
                        )
                        output_file = re.search("(.*?)" + genome_pattern, genome_file)
                        for output_lines in output_lines_lst:
                            with open(
                                output_file.group(1)
                                + "_"
                                + output_lines[1]
                                + "_extraction_RAW.fas",
                                "w",
                            ) as file:
                                first = True
                                for line in output_lines[0]:
                                    file.write(">")
                                    if first:
                                        first = False
                                        file.write(line[0])
                                        file.write("\n")
                                        file.write(line[1])
                                        file.write("\n")
                                    else:
                                        for inner_line in line:
                                            file.write(inner_line[0])
                                            file.write("\n")
                                            file.write(inner_line[1])
                                            file.write("\n")
        if not found and (ID.group() is not None):
            logger.warning(
                f"Not able to find any genome file with '{ID.group(1)}'. Genome and Blast unique output file must have a 'Genus_species_assemblyID' structure or the name provided by NCBI when downloaded."
            )


@click.command()
@click.option(
    "--genomic",
    help="Extract unique hits from genomic data stored at Genomic folder or specify sequence type if using --directory argument.",
    type=click.STRING,
    is_flag=True,
    flag_value="genomic",
)
@click.option(
    "--rna",
    help="Extract unique hits from rna data stored at Rna folder or specify sequence type if using --directory argument.",
    type=click.STRING,
    is_flag=True,
    flag_value="rna",
)
@click.option(
    "--protein",
    help="Extract unique hits from protein data stored at Protein folder or specify sequence type if using --directory argument.",
    type=click.STRING,
    is_flag=True,
    flag_value="protein",
)
@click.option("--directory", help="Path to custom folder", type=click.STRING)
@click.option(
    "--blast_pattern",
    help="Custom pattern to find BLAST unique output files.",
    default="unique.tsv",
    type=click.STRING,
    show_default=True,
)
@click.option(
    "--genome_pattern",
    help="Pattern to find genome files.",
    default=".f[a,n]a",
    type=click.STRING,
    show_default=True,
)
@click.option(
    "--in_len",
    help="Number of sites extracted upstream and downstream from the blast hit.",
    default=10000,
    type=click.INT,
    show_default=True,
)
@click.option(
    "--query_seqs",
    help="Use it when you want to attach some query sequences according to BLAST results for future alignments",
    is_flag=True,
)
@click.option(
    "--query_seqs_num",
    help="Maximum number of query sequences extracted. Use it when you want to attach some sequences to genomic sequences for future alignments.",
    default=5,
    type=click.INT,
    show_default=True,
)
@click.option(
    "--email",
    help="Mandatory when you want to download some sequences to complete the RAW files for future alignments. As in get_query_seqs.",
    type=click.STRING,
)
def main(
    genomic,
    rna,
    protein,
    directory,
    blast_pattern,
    genome_pattern,
    in_len,
    query_seqs,
    query_seqs_num,
    email,
):
    """Extract RAW sequences from whole genome derived data. Implemented for Genus_species_assemblyID or GC[A,F]_ID. Both Blast output and genome data files must have the same structure."""
    path = os.getcwd()
    data_type_lst = [genomic, rna, protein]

    if directory is not None:
        found = False
        for data_type in data_type_lst:
            if data_type is not None:
                found = True
                logger.info(" Running directory and %s argument.\n" % data_type)
                generate_output_FASTA_file(
                    directory + "*",
                    data_type,
                    blast_pattern,
                    genome_pattern,
                    query_seqs,
                    query_seqs_num,
                    in_len,
                    email,
                )
        if not found:
            logger.error(". You must use genomic, rna or protein arguments.")
    else:
        found = False
        for data_type in data_type_lst:
            if data_type is not None:
                found = True
                logger.info(f"Running {data_type} argument.")
                generate_output_FASTA_file(
                    path + "/../Data/" + data_type.capitalize() + "/*",
                    data_type,
                    blast_pattern,
                    genome_pattern,
                    query_seqs,
                    query_seqs_num,
                    in_len,
                    email,
                )
        if not found:
            logger.error(f"You must use genomic, rna, protein arguments.")
