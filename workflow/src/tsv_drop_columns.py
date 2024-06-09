import argparse
import pandas as pd
import sys

# parse command line arguments
parser = argparse.ArgumentParser(description = "Removes given columns from tsv-file.")
parser.add_argument("--input", type = str, help = "input file (tsv)")
parser.add_argument("--columns", type = str, help = "column numbers to be deleted. Multiple column numbers separated by comma, ranges given by hyphen, e.g. 4,7,10-12")
parser.add_argument("--output", type = str, help = "output file (tsv)")
args = parser.parse_args()

if args.input is None or args.columns is None or args.output is None:
    parser.print_help()
    sys.exit()

### xlsx to tsv conversion
data = pd.read_csv(args.input, header = None, index_col = None, sep = '\t')

column_indices = []
for column_range in args.columns.split(','):
	if '-' in column_range:
		column_range_positions = column_range.split('-')
		for column_index in range(int(column_range_positions[0]) - 1, int(column_range_positions[1])):
			column_indices.append(column_index)
	else:
		column_indices.append(int(column_range) - 1)

data.drop(columns = data.columns[column_indices], inplace = True)

data.to_csv(args.output, encoding = 'utf-8', header = False, index = False, sep = '\t')
