# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Filtering of genotype data against missing alleles/uniallelic markers

import argparse
import sys
import numpy as np
import pandas as pd

# parse command line arguments
parser = argparse.ArgumentParser(description = 'Filters SNP genotype data.')
parser.add_argument('--genotype', type = str, help = 'Genotype data (tsv file)')
parser.add_argument('--header_columns', type = int, help = 'number of non-genotype columns')
parser.add_argument('--output', type = str, help = 'output file (tsv)')
args = parser.parse_args()

if args.genotype is None or args.output is None or args.header_columns is None:
    parser.print_help()
    print("\nAn input for the arguments: genotype, header_columns and output is required.")
    sys.exit()

# read genotyping data
genotype_data = pd.read_csv(args.genotype, header = None, index_col = None, delimiter = '\t', low_memory=False)

total_marker_count = len(genotype_data.index[args.header_columns:])

print('Number of markers: \t' + str(total_marker_count))

alternate_allele_column_index = [index for index, element in enumerate(genotype_data.loc[0].tolist()) if 'Alternate allele' in element][0]

filtered_marker_rows = []

for marker_index in range(1, len(genotype_data[0])):
    # filter missing alternative alleles
    if '.' in genotype_data.iloc[marker_index, alternate_allele_column_index]:
        filtered_marker_rows.append(marker_index)
        continue
    
    # filter multiple alternative alleles
    if ';' in genotype_data.iloc[marker_index, alternate_allele_column_index]:
        filtered_marker_rows.append(marker_index)
        continue
    
    marker_allele_counts = pd.value_counts(genotype_data.iloc[marker_index, args.header_columns:])
    
    # filter marker with only one known allele
    if len(marker_allele_counts) <= 1 or (('.' in marker_allele_counts.index and len(marker_allele_counts) == 2)):
        filtered_marker_rows.append(marker_index)
        continue

genotype_data.drop(labels = filtered_marker_rows, axis = 0, inplace = True)

print('Number of markers after filtering: ' + str(len(genotype_data.index[args.header_columns:])))
print('Number of genotypes after filtering: ' + str(len(genotype_data.iloc[0])))

# output filtered genotype data
genotype_data.to_csv(args.output, encoding = 'utf-8', header = False, index = False, sep = '\t')
