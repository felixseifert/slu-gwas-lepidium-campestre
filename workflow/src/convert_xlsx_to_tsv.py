# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Extract Excel sheet to tab-separated file (tsv)

import argparse
import pandas as pd
import sys

# parse command line arguments
parser = argparse.ArgumentParser(description = "Export sheet from xlsx-file into tsv-file.")
parser.add_argument("--xlsx", type = str, help = "Excel file (xlsx)")
parser.add_argument("--sheet", type = str, help = "xlsx sheet title")
parser.add_argument("--output", type = str, help = "output file (tsv)")
args = parser.parse_args()

if args.xlsx is None or args.sheet is None or args.output is None:
    parser.print_help()
    sys.exit()

### xlsx to tsv conversion
data_xls = pd.read_excel(args.xlsx, args.sheet, header = None, index_col = None)
data_xls.to_csv(args.output, encoding = 'utf-8', header = False, index = False, sep = '\t')
