import argparse
import sys

import pandas as pd

# parse command line arguments
parser = argparse.ArgumentParser(description='Reformats clumpp result file for pophelper')
parser.add_argument('--input', type=str, help='clumpp result file')
parser.add_argument('--output', type=str, help='output file')
args = parser.parse_args()

if args.input is None or args.output is None:
    parser.print_help()
    sys.exit()

clumpp_df = pd.read_csv(args.input, header=None, index_col=0, sep=" ", skipinitialspace=True)
cluster_order = clumpp_df.sum()[4:len(clumpp_df.columns)-1].sort_values(ascending=False).index

clumpp_df[1] = clumpp_df[1].astype(str) + ':'
clumpp_df[2] = 1
clumpp_reformat_df = clumpp_df.reindex([1] + cluster_order.values.tolist() + [2], axis='columns')

clumpp_reformat_df.to_csv(args.output, encoding = 'utf-8', header=None, index=None, sep=' ')