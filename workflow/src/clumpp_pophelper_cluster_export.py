# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Extraction of cluster index with highest supprot from pophelper-formatted CLUMPP result

import argparse
import os
import sys

import pandas as pd

# parse command line arguments
parser = argparse.ArgumentParser(description='Extracts most probable cluster index pophelper-formatted clumpp result')
parser.add_argument('--cluster', type=str, help='pophelper-formatted clumpp input file')
parser.add_argument('--label', type=str, help='sample labels')
parser.add_argument('--output', type=str, help='output file')
args = parser.parse_args()

if args.cluster is None or args.output is None:
    parser.print_help()
    sys.exit()

cluster_df = pd.read_csv(args.cluster, sep=' ', header=None, index_col=None)
label_df = pd.read_csv(args.label, sep=' ', header=None, index_col=None)
 
result_df = cluster_df[cluster_df.columns[1:-1]].idxmax(axis=1)
result_df.index = label_df[0].values

result_df.to_csv(args.output, header=None, index=False, sep='\t')