# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Generation of clustered heatmap from Plink additive recoding

import argparse
import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# parse command line arguments
parser = argparse.ArgumentParser(description = 'Generates clustered heatmap from plink additive recoding (--recodeA).')
parser.add_argument('--genotype', type = str, help = 'plink Genotype dataset(s) (tsv file)', nargs = '*')
parser.add_argument('--output', type = str, help = 'output file prefix (without ending)')
args = parser.parse_args()

if args.genotype is None or args.output is None:
    parser.print_help()
    sys.exit()

first_file_flag = True
for plink_recodea_file in args.genotype:
    if first_file_flag:
        plink_recodea_df = pd.read_csv(plink_recodea_file, delimiter=" ", header=0, index_col=0)
        first_file_flag = False
    else:
        plink_recodea_additional_df = pd.read_csv(plink_recodea_file, delimiter=" ", header=0, index_col=0)
        plink_recodea_df = pd.concat([plink_recodea_df, plink_recodea_additional_df])

data_df = plink_recodea_df[plink_recodea_df.columns[5:]]
data_df = data_df.fillna(-3)

clustermap = sns.clustermap(data_df, cmap=sns.diverging_palette(20, 220, as_cmap=True), xticklabels=False, yticklabels=False)
clustermap.cax.set_visible(False)
ax = clustermap.ax_heatmap
ax.set_xlabel('Marker alleles')
ax.set_ylabel('Lines')

clustermap.savefig(args.output + '.png', dpi=300)
#clustermap.savefig(args.output + '_mds_matrix_plot.svg', dpi=300)