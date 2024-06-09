# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Removal of sample ids from dataset

import argparse
import sys
import numpy as np
import pandas as pd

# parse command line arguments
parser = argparse.ArgumentParser(description = 'removes sample ids from dataset')
parser.add_argument('--input', type = str, help = 'input fasta file')
parser.add_argument('--samples', type=str, help = 'sample ids to be removed tsv file')
parser.add_argument('--output', type = str, help = 'output fasta file')
args = parser.parse_args()

if args.input is None or args.output is None or args.samples is None:
    parser.print_help()
    print("\nAn input for the arguments: input and output as well as samples is required.")
    sys.exit()

# load sample to be removed into dict
samples_remove = list()
header_flag = True
with open(args.samples) as samples_file:
    for line in samples_file:
        samples_remove.append(line)

# export filtered sequence ids and save to output
outfile = open(args.output, 'w')

with open(args.input, 'r') as infile:
    for line in infile:
        dataset_parts = line.split('\t', 1)
        sample_id = str(dataset_parts[0]).replace("Lep", "LEP")
        
        if not sample_id in samples_remove:
            outfile.write(sample_id + '\t' + dataset_parts[1])

outfile.close()