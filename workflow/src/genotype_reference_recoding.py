# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Geneotype data reference recoding based on reference

import argparse
import sys
import pandas as pd
 
# parse command line arguments
parser = argparse.ArgumentParser(description = 'Replace position on scaffold with assembly reference position')
parser.add_argument('--genotype', type = str, help = 'genotype data input file (tsv file)')
parser.add_argument('--agp', type = str, help = 'assembly genome path input file (tsv file)')
parser.add_argument('--output', type = str, help = 'output file (csv file)')
args = parser.parse_args()

if args.genotype is None or args.agp is None or args.output is None:
    parser.print_help()
    sys.exit()

# read tsv data
genotype_df = pd.read_csv(args.genotype, delimiter='\t', header=None, index_col=0).T

agp_df = pd.read_csv(args.agp, delimiter = '\t', header=None, index_col=None)
agp_df.columns = ['chromosome', 'chromosome_start', 'chromosome_end', 'chromosome_rank', 'type', 'sequence', 'sequence_start', 'sequence_end', 'orientation']
agp_df.index = agp_df['sequence']

# remove gaps
agp_df = agp_df[agp_df['sequence'] != '100']


# replace anchored sequences with chromosome and position
for index in agp_df.index:
    chromosome = agp_df.loc[index, 'chromosome']
    chromosome_start = int(agp_df.loc[index, 'chromosome_start'])
    sequence = agp_df.loc[index, 'sequence']
    sequence_end = int(agp_df.loc[index, 'sequence_end'])
    orientation = agp_df.loc[index, 'orientation']
    
    for genotype_index in genotype_df[genotype_df['Reference'] == sequence].index.values.tolist():
        marker_sequence_position = int(genotype_df.loc[genotype_index, 'Position'])
        
        marker_reference_position = (chromosome_start + marker_sequence_position)
        if orientation == '-':
            marker_reference_position = (chromosome_start + sequence_end - marker_sequence_position)
        
        genotype_df.loc[genotype_index, 'Position'] = marker_reference_position
        genotype_df.loc[genotype_index, 'Reference'] = chromosome

# output transposed data
genotype_df.to_csv(args.output, encoding = 'utf-8', header=True, index=False, sep='\t')
