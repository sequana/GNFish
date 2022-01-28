#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import urllib
from Bio import Entrez
from Bio import GenBank
from Bio import SeqIO
import re
import argparse
import csv

 
##Arguments
parser = argparse.ArgumentParser(description='Downloads a set of protein for query for BLAST searches. Alternatively you can dowload nucleotides sequences using IDs')
parser.add_argument('email', help='sets mandatory e-mail for NCBI searches', type=str)
parser.add_argument('query', help='sets file path with the query', type=argparse.FileType('r'))
parser.add_argument('--nucleotide', help='for nucleotide downloading', action='store_true')
parser.add_argument('--refine', help='adds filter or field tags to the query. Constant value refseq[filter]. Follow constant value structure for your custom filters, begin with "AND".', nargs='?', const=' AND refseq[filter]', type=str, default='')
parser.add_argument('--curated', help='only donwloads sequences with protein name included in Protein Feature. Just for protein searches', action='store_true')
parser.add_argument('--retmax', help='sets number of searched NCBI records. Default value equal to 200',nargs='?', const=200, type=int, default=200)
args = parser.parse_args()
query_file = args.query
email = args.email
nucleotide = args.nucleotide
retmax = args.retmax
refine = args.refine
curated = args.curated
db='protein'


if nucleotide:
    db='nucleotide'
path = os.getcwd()
query_lst = query_file.read().split('\n')

def create_directory(directory):
    directory = path +'/'+ directory
    if not os.path.isdir(directory):
        os.mkdir(directory)
    else:
        print('WARNING. %s directory already exits. New data will be added to previous one.\n' %directory)

def print_running_info(query, retmax):
    print('INFO. Queries must follow gene name followed by filters between parentheses. If no filters type () after gene name. For nucletoide add --refine biomol_mrna[PROP] for proper transcrit download.\n')
    print('INFO. Running %s query, %s data type, and retmax = %s. If using --curated note that the final number of sequences can be less than retmax.\n' %(query, db, retmax))

def get_assembly_summary(id_num):
    esummary_handle = Entrez.esummary(db=db, id=id_num, report="full")
    esummary_record = Entrez.read(esummary_handle, validate=False)
    return esummary_record

def check_protein_in_all_values(features, protein):
    found = False
    for key, value in features.items():
        if type(value) is dict:
            check_protein_in_all_values(features, protein)
        # elif type(value) is list:
        #     for i in range(len(value)):
        #          get_all_values(value[i])
        else:
            if re.search('product', str(value)):
                value = value[0]['GBQualifier_value']
                if re.search(protein, value):
                    found = True
            if not found:
                if re.search('note', str(value)):
                    value = value[0]['GBQualifier_value']
                    if re.search(protein, value):
                        found = True
    return found

def get_seqs_data(term, retmax, name):
    Entrez.email = email
    handle = Entrez.esearch(db=db, term=term, retmax=retmax, sort='Significance')
    record = Entrez.read(handle)
    handle.close()
    rows = []
    try:
        error_sentence=record['ErrorList']['PhraseNotFound']
        print('ERROR. %s from %s not found and search result can be unpredictable. Check %s on "https://www.ncbi.nlm.nih.gov/" protein database or correct %s term.\n' %(error_sentence, term, term, error_sentence))
    except KeyError:
        ids = record['IdList']
        for id_num in ids:
            flag=True
            if curated and db=='protein':
                flag=False
                handle= Entrez.efetch(db=db, id=id_num, retmode="xml")
                parser = Entrez.read(handle, validate=True, escape=False)
                features = parser[0]['GBSeq_feature-table']
                for i in range(len(features)):
                    found = check_protein_in_all_values(features[i], name)
                    if found:
                        flag=True
            if flag:
                handle = Entrez.efetch(db=db, id=id_num,rettype='fasta', retmode="text")
                record = handle.read()
                rows.append(record)
    with open(path +'/'+ name+'_query_seqs.fas', 'w') as file:
        for row in rows:
            file.write(row)
            

path = path +'/../Data'
create_directory('Query_seqs')
path = path +'./Query_seqs'
for query in query_lst:
    if query=='':
        continue
    else:
            term=query
            if refine:
                term = term+refine
            a = re.search('(.*)?\(', term)
            if a:
                name = re.sub(' ', '_',a.group(1))
                print_running_info(term, retmax)
                get_seqs_data(term, retmax, name)
            else:
                print('WARNING. %s query is not correct. Must follow Gene name followed by filters between parentheses. If no filters type () after protein name. Skipt to next query.' %(query))



