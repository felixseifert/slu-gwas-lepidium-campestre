# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Merge structure result files for CLUMPP analysis

import argparse
import os
import re
import sys

# parse command line arguments
parser = argparse.ArgumentParser(description='Transposes tab-separated file')
parser.add_argument('--structure', type=str, help='directory containing structure result files')
parser.add_argument('--prefix', type=str, help='shared prefix of structure result files')
parser.add_argument('--output_ind', type=str, help='clumpp indfile')
parser.add_argument('--output_pop', type=str, help='clumpp popfile')
args = parser.parse_args()

if args.prefix is None or args.output_ind is None or args.output_pop is None:
    parser.print_help()
    sys.exit()

structure_directory = args.structure
if structure_directory == None:
    structure_directory = '.'

individual_number_data_pattern = re.compile('([0-9]*)\sindividuals')
structure_individual_data_pattern = re.compile('(\s*[0-9]*)\s*(\([0-9]*\))\s*:\s(.*)')
structure_population_data_pattern = re.compile('([01]\.[0-9]{3}\s\s){2,}')

try:
    output_ind_fp = open(args.output_ind, 'w')
    output_pop_fp = open(args.output_pop, 'w')

    number_individuals = "1"
    
    for file in os.listdir(structure_directory):
        if file.startswith(args.prefix):
            with open(structure_directory + "/" + file) as current_structure_replicate_file:
                cluster_information_flag = False

                for line in current_structure_replicate_file:
                    individual_number_match = individual_number_data_pattern.search(line)
                    individual_match = structure_individual_data_pattern.search(line)
                    population_match = structure_population_data_pattern.match(line)
                    
                    if individual_number_match:
                        number_individuals = individual_number_match.group(1)
                    elif individual_match:
                        line_id = individual_match.group(1).strip()
                        missing = individual_match.group(2).strip()
                        inferred_clusters = individual_match.group(3).lstrip()
                        
                        output_ind_fp.write(line_id + ' ' + line_id + ' ' + missing + ' 1 : ' + inferred_clusters + '\n')
                    elif population_match:
                        output_pop_fp.write('1: ' + population_match.group(0) + str(number_individuals))
            
            output_ind_fp.write('\n')
            output_pop_fp.write('\n')

except Exception as error_message:
    print('An error occured: ' + str(error_message))
finally:
    output_ind_fp.close()
    output_pop_fp.close()