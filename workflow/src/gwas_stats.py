# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Summary generation of GWAS analyses with different PCA covariate number

import sys
import pandas as pd

args = sys.argv[1:]

if len(args) == 0:
    parser.print_help()
    print("\nOne or more input files needs to be specified.")
    sys.exit()


print("trait\tcovar_type\tcovar_number\tmin_pvalue\tcount_significant\tnumber_marker")

for filename in args:
    data = pd.read_excel(filename)

    count_significant = (data["FDR_BY"] < 0.05).sum()
    min_pvalue = min(data["FDR_BY"])
    number_marker = len(data)

    filename_parts = filename.replace("\\", "/").split("/")[-1].split(".")[0].split("_")

    print(f"{' '.join(filename_parts[0:-2])}\t{filename_parts[-2]}\t{filename_parts[-1]}\t{min_pvalue}\t{count_significant}\t{number_marker}")
