import argparse
import sys
import pandas as pd

# parse command line arguments
parser = argparse.ArgumentParser(description = 'Transposes tab-separated file')
parser.add_argument('--input', type = str, help = 'input file (tsv file)')
parser.add_argument('--output', type = str, help = 'output file (csv file)')
args = parser.parse_args()

if args.input is None or args.output is None:
    parser.print_help()
    sys.exit()

# read tsv data
data = pd.read_csv(args.input, delimiter = '\t', header = None, index_col = None)

# output transposed data
data.transpose().to_csv(args.output, encoding = 'utf-8', header = False, index = False, sep = '\t')
