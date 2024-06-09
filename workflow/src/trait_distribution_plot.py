#!/bin/env python3

import argparse
import sys

import pandas as pd
import seaborn as sns

# parse command line arguments
parser = argparse.ArgumentParser(description = "Generates trait distribution plot.")
parser.add_argument("--input", type = str, help = "input file (tsv)")
parser.add_argument("--trait_label", type = str, help = "trait axis label")
parser.add_argument("--output", type = str, help = "output filename")
args = parser.parse_args()

if args.input is None or args.trait_label is None or args.output is None:
    parser.print_help()
    sys.exit()


trait_df = pd.read_csv(args.input, sep="\t")
trait_column = trait_df.columns[2]

trait_plot = sns.displot(trait_df[trait_column], kind="kde")
trait_plot.set_axis_labels(args.trait_label.replace("_", " "))
trait_plot.figure.savefig(args.output)
