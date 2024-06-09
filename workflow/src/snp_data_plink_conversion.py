# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Conversion of SNP data to Plink format

import argparse
import sys
import numpy as np
import pandas as pd

# parse command line arguments
parser = argparse.ArgumentParser(description = 'Converts SNP genotyping data for Plink.')
parser.add_argument('--genotype', type = str, help = 'Genotype data (tsv file)')
parser.add_argument('--output', type = str, help = 'output file prefix without file-endings')
args = parser.parse_args()

if args.genotype is None or args.output is None:
    parser.print_help()
    print("\nAn input for the arguments: genotype and output is required.")
    sys.exit()

# read genotyping data
genotype_data = pd.read_csv(args.genotype, header = 0, index_col = None, delimiter = '\t', low_memory=False)

fam_output = open(args.output + '.fam', 'w')
map_output = open(args.output + '.map', 'w')
ped_output = open(args.output + '.ped', 'w')

for genotype_index in range(6, len(genotype_data.iloc[0])):
    fam_output.write('0\t' + genotype_data.columns[genotype_index] + '\t0\t0\t0\t-9\n')
    ped_output.write('0\t' + genotype_data.columns[genotype_index] + '\t0\t0\t0\t-9')
    
    for marker_index in range(0, len(genotype_data.index)):
        alleles = '0\t0'
        
        if(len(genotype_data.iloc[marker_index, genotype_index]) > 1):
            alleles = genotype_data.iloc[marker_index, genotype_index].replace('/', '\t').replace('.', '0')
        
        ped_output.write('\t' + alleles)
    
    ped_output.write('\n')

for marker_index in range(0, len(genotype_data.index)):
    map_output.write(str(genotype_data.iloc[marker_index, 0]) + '\t' + str(genotype_data.iloc[marker_index, 2]) + '\t0\t' + str(genotype_data.iloc[marker_index, 1]) + '\n')

fam_output.close()
map_output.close()
ped_output.close()
