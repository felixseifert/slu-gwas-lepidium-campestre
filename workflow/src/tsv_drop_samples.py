import argparse
import pandas as pd
import sys

# parse command line arguments
parser = argparse.ArgumentParser(description = "Removes given samples from tsv-file.")
parser.add_argument("--input", type = str, help = "input file (tsv)")
parser.add_argument("--samples", type = str, help = "name of file containing sample names (one per row)")
parser.add_argument("--output", type = str, help = "output file (tsv)")
args = parser.parse_args()

if args.input is None or args.samples is None or args.output is None:
    parser.print_help()
    sys.exit()

### xlsx to tsv conversion
data = pd.read_csv(args.input, header = None, index_col = None, sep = '\t')

samples = pd.read_csv(args.samples, header=None,  index_col = None)
sample_names = samples[0].tolist()

data_filtered = data[~data[0].isin(sample_names)]

data_filtered.to_csv(args.output, encoding = 'utf-8', header = False, index = False, sep = '\t')
