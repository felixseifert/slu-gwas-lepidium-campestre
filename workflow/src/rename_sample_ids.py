# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Renaming of sample ids

import argparse
import sys
import numpy as np
import pandas as pd

# parse command line arguments
parser = argparse.ArgumentParser(description = 'replaces assembly sequence id by scaffold/contig name')
parser.add_argument('--input', type = str, help = 'input fasta file')
parser.add_argument('--translation', type=str, help = 'sample id translation tsv file')
parser.add_argument('--output', type = str, help = 'output fasta file')
args = parser.parse_args()

if args.input is None or args.output is None or args.translation is None:
    parser.print_help()
    print("\nAn input for the arguments: input and output as well as sample translation is required.")
    sys.exit()

# load sample translation into dict
sample_translation = {}
header_flag = True
with open(args.translation) as translation_file:
    for line in translation_file:
        if header_flag:
            # skip header line
            header_flag = False
        else:
            (key, value) = line.split()
            sample_translation[key] = value

# translate sequence ids and save to output
outfile = open(args.output, 'w')

with open(args.input, 'r') as infile:
    for line in infile:
        dataset_parts = line.split('\t', 1)
        sample_id = str(dataset_parts[0]).replace("Lep", "LEP")
        
        if sample_id in sample_translation:
            outfile.write(sample_translation[sample_id] + '\t' + dataset_parts[1])
        else:
            outfile.write(line)

outfile.close()