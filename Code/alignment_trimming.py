#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 21:41:03 2021

@author: hectorlorente
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 16:32:14 2021

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
parser = argparse.ArgumentParser(description='Alignment trimming using trimAl.')
parser.add_argument('trimal_path', help='path to trimAl directory', type=str)
parser.add_argument('input_file', help='alignment (FASTA file with prot or nucl aligned seqs)', type=str)
parser.add_argument('output_file', help='output file name.', nargs='?', type=str)
parser.add_argument('--trimal_command', help='command to run trimal options. Default -gt 0.1', nargs='?', const='-gt 0.1', type=str, default='-gt 0.1')
args = parser.parse_args()
trimal_path = re.sub('\'', '', args.trimal_path)
input_file = re.sub('\'', '', args.input_file)
output_file = re.sub('\'', '', args.output_file)
trimal_command = args.trimal_command
output_suffix = args.output_suffix
path = os.getcwd()


trimalCommand = [trimal_path + 'trimal -s ', input_file, '-out', output_file, trimal_command]
subprocess.run(trimalCommand)

