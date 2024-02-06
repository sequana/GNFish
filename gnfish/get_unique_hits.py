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

##Arguments
parser = argparse.ArgumentParser(description='Gets unique hits from BLAST output files based on genomes IDs.')
parser.add_argument('--genomic', help='extracts unique hits from genomic data stored at Genomic folder. Use directory for custom folder.', nargs='?', const='genomic', type=str)
parser.add_argument('--rna', help='extracts unique hits from rna data stored at Rna folder. Use directory for custom folder.', nargs='?', const='rna', type=str)
parser.add_argument('--protein', help='extracts unique hits from protein data stored at Protein folder. Use directory for custom folder.', nargs='?', const='protein', type=str)
parser.add_argument('--directory', help='sets path to custom folder',nargs='?', type=str)
parser.add_argument('--pattern', help='custom pattern to find Blast output files. Default ".tsv".', nargs='?', default='.tsv', const='.tsv', type=str)
args = parser.parse_args()
genomic = args.genomic
rna = args.rna
protein = args.protein
directory = args.directory
pattern = args.pattern
get_unique = args.get_unique
path = os.getcwd()
data_type_lst=[genomic, rna, protein]

def open_TSV_file(input_file):
    rows = dict()
    with open(input_file) as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        return[row for row in tsv_reader]
    
def select_files(path):
    return list_files.list_files_method(list_files,path)

def get_files(path):
    files = []
    found = False
    for folder in select_files(path):
        try:
            for i in range(len(folder[1])):
                file = re.search('(.*)' + pattern, folder[1][i])
                if file and (not re.search('_unique.tsv', file.group())) and (not re.search('_filtered.tsv', file.group())):
                    in_file = re.sub('\\\\', '/', file.group())
                    files.append(in_file)
                    found = True
                    # print(file.group())
        except IndexError:
            continue
    if not found:
        print('ERROR. No file. Try another directory or pattern. Current pattern extension is "%s".'%pattern)
        sys.exit()
    else:
        return files

def get_unique_hits(rows):
    new_rows = []
    rows_1 =[]
    first = True
    for row in rows:
        if first:
            new_rows.append(row)
            rows_1.append(row[1])
            first = False
        else:
            if not (row[1] in rows_1):
                new_rows.append(row)
                rows_1.append(row[1])
    return new_rows
        
def generate_output_file(path):
    tsv_file_lst = get_files(path)
    for tsv_file in tsv_file_lst:
        new_rows = []
        exten = ''
        print('Extracting hits from %s.\n' %tsv_file)
        rows = open_TSV_file(tsv_file)
        # if len_treshold:
        #     print('Filtering hits from %s.\n' %tsv_file[0])
        #     exten = '_filtered'
        #     if get_unique:
        #         new_rows = get_unique_hits(rows)
        #         exten = exten + '_unique'
        # elif get_unique:
        new_rows = get_unique_hits(rows)
        exten = 'unique.tsv'
        # else:
        #     print('ERROR. You must use "get_unique" or "filter" arguments.')
        #     sys.exit()
        # print('aqui')
        # print(re.sub(pattern, '_'+ exten, tsv_file))
        with open(re.sub(pattern, '_'+ exten, tsv_file), 'w') as file:
            for row in new_rows:
                for i in range(len(row)):
                    file.write(row[i] + '\t')
                file.write('\n')
                
if directory is not None:
    print('INFO. Running directory argument. Blast files must have .tsv extension.\n')
    generate_output_file(directory + '*')
else:
    found = False
    for data_type in data_type_lst:
        if data_type is not None:
            found = True
            print('INFO. Running %s argument. Blast files must have %s extension.\n' %(data_type, pattern))
            generate_output_file(path+'/../Data/'+ data_type.capitalize()+'/*')
    if not found:
        print('ERROR. You must use genomic, rna, protein or directory arguments.')


