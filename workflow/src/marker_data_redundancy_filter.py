# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Filtering against redundant markers and genotypes/markers with high missingness

import argparse
import itertools
import math
import os
import sys

import numpy as np
import pandas as pd

max_genotype_identity = 0.99
max_marker_identity = 0.99
max_missing_genotype_fraction = 0.75
max_missing_marker_fraction = 0.1

# parse command line arguments
parser = argparse.ArgumentParser(description = 'Filters plink dataset for redundant genotypes/SNPs.')
parser.add_argument('--input', type = str, help = 'Genotype dataset (tsv file)')
parser.add_argument('--genotype', type = float, help = 'maximum genotype identity (default: ' + str(max_genotype_identity))
parser.add_argument('--marker', type = float, help = 'maximum marker identity (default: ' + str(max_marker_identity))
parser.add_argument('--missing_genotype', type = float, help = 'maximum missing genotype fraction (default: ' + str(max_missing_genotype_fraction))
parser.add_argument('--missing_marker', type = float, help = 'maximum missing marker fraction (default: ' + str(max_missing_marker_fraction))
parser.add_argument('--output', type = str, help = 'output file prefix (without ending)')
args = parser.parse_args()

if args.input is None or args.output is None:
    parser.print_help()
    sys.exit()

if args.genotype is not None:
    max_genotype_identity = args.genotype

if args.marker is not None:
    max_marker_identity = args.marker

if args.missing_genotype is not None:
    max_missing_genotype_fraction = args.missing

if args.missing_marker is not None:
    max_missing_marker_fraction = args.missing

# read genotyping data
genotype_dataset_df = pd.read_csv(args.input, header = 0, index_col = 1, delimiter = ' ')
genotype_data_df = genotype_dataset_df[genotype_dataset_df.columns[5:]].copy()
    
# filter missing marker
filtered_markers = []
for marker in genotype_data_df.columns:
    if genotype_data_df[marker].isna().sum() > (max_missing_marker_fraction * len(genotype_data_df[marker])):
        filtered_markers.append(marker)

genotype_dataset_df.drop(labels = filtered_markers, axis = 1, inplace = True)
genotype_data_df.drop(labels = filtered_markers, axis = 1, inplace = True)

# filter missing genotype
filtered_genotypes = []
for genotype in genotype_data_df.index:
    if genotype_data_df.loc[genotype].isna().sum() > (max_missing_genotype_fraction * len(genotype_data_df)):
        filtered_genotypes.append(genotype)

genotype_dataset_df.drop(labels = filtered_genotypes, axis = 0, inplace = True)
genotype_data_df.drop(labels = filtered_genotypes, axis = 0, inplace = True)

if not os.path.exists(args.output):
    os.mkdir(args.output)

marker_count_raw = len(genotype_data_df)
genotype_count_raw = len(genotype_data_df.columns)

# remove marker redundancy
marker_correlation_matrix = np.abs(genotype_data_df.corr())
correlated_marker_features = []

file_object_marker_redundancy = open(args.output + '_redundant_marker.csv', 'w')

for i in range(len(genotype_data_df.columns) - 2):
    for j in range(i):
        if abs(marker_correlation_matrix.iloc[i, j]) > max_marker_identity:
            if marker_correlation_matrix.columns[i] not in correlated_marker_features:
                file_object_marker_redundancy.write(marker_correlation_matrix.columns[i] + '\t' + marker_correlation_matrix.columns[j] + '\n')
                correlated_marker_features.append(marker_correlation_matrix.columns[i])

file_object_marker_redundancy.close()

genotype_dataset_df.drop(labels = correlated_marker_features, axis = 1, inplace = True)
genotype_data_df.drop(labels = correlated_marker_features, axis = 1, inplace = True)

# remove genotype redundancy
genotype_data_df_transposed = genotype_data_df.T
genotype_correlation_matrix = np.abs(genotype_data_df_transposed.corr())
correlated_genotype_features = []

file_object_genotype_redundancy = open(args.output + '_redundant_genotype.csv', 'w')

for i in range(len(genotype_data_df_transposed.columns) - 2):
    for j in range(i):
        if abs(genotype_correlation_matrix.iloc[i, j]) >= max_genotype_identity:
            if genotype_correlation_matrix.columns[i] not in correlated_genotype_features:
                file_object_genotype_redundancy.write(genotype_correlation_matrix.columns[i] + '\t' + genotype_correlation_matrix.columns[j] + '\n')
                correlated_genotype_features.append(genotype_correlation_matrix.columns[i])

file_object_genotype_redundancy.close()

genotype_dataset_df.drop(labels = correlated_genotype_features, axis = 0, inplace = True)

genotype_dataset_df.to_csv(args.output + '_non_redundant.raw', sep = ' ', encoding = 'utf-8')
