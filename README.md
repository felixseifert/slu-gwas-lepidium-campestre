
# SLU Lepidium campestre GWAS

This repository contains the Snakemake pipeline along with all custom code to run the SLU Lepidium campestre GWAS.
The pipeline is implemented in Snakemake and uses custom Python and R scripts and bash scripts.

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.81-blue.svg)](https://snakemake.github.io)
[![Conda](https://img.shields.io/badge/Conda-brightgreen.svg?style=flat-square)](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
![License](https://img.shields.io/badge/license-MIT-black.svg)

Note: This repository only covers the code of the pipeline.

## Installation

### Installing via Conda 

In order to use this workflow, you need to clone this repository and have Conda [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed (to install `conda`, please follow the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)).

We will set a Conda environment for running Snakemake which will be created as follows:

```bash
# Create a new conda environment with Snakemake
# and its dependencies
conda create -c bioconda -n snakemake -f workflow/envs/snakemake.yaml

# Activate the Snakemake conda environment
conda activate snakemake
```

### Running the Snakemake pipeline

when the Snakemake conda environment is activated you execute the pipeline as follows:

```bash
# Run the GWAS pipeline (adapt the number of cores to your machine)
snakemake --cores 4 --use-conda -s workflow/Snakemake
```

Upon execution the following directories will be generated within the location you execute the Snakemake pipeline:

```
path/to/your/work_dir/
├── analysis/
│   ├── data/       Reformatted input data
│   ├── genotype/   Genotype data of the various processing steps
│   ├── gwas/       GWAS analysis results for different filterings of the dataset and per phenotype
│   ├── marker_anchoring/   Results from anchoring of SNP markers to the reference genome
│   ├── phenotype/  Reformatted and adjusted phenotype data
│   ├── plots/      Plots generated for the GWAS analysis
│       ├── gwas/       Manhattan and QQ-plots  of different filterings of the datasets and per phenotype
│       ├── heatmaps/   Heatmaps of the genotype information at different filtering steps
│       ├── pca/        PCA plots and variants applying structure populations of various k
│       ├── phenotype/  Distribution plots of adjusted phenotype values
│       ├── population_structure/   Population structure analysis plots of various k
│   ├── population_structure/   Analysis results for population structure analysis (PCA, structure)
├── benchmark/        benchmark-files created for the rules of the Snakemake pipeline
├── logs/             Log-files created for the rules of the Snakemake pipeliine
```

### Directory Structure 
```
path/to/your/work_dir/
├── workflow/
│   ├── envs/         directory containing Conda environemnt definitions (yaml files) for the pipeline 
│   ├── src/          directory containing custom scripts (Python and R) for the pipeline 
│   ├── Snakefile     The Snakemake pipeline for the GWAS analysis
├── LICENSE           The MIT license for this project
├── README.md         This README fil
```

## Authors

The GWAS code was developed by Felix Seifert, cropSeq bioinformatics.

## License
The pipeline and custom code is licensed under the [MIT](LICENSE) license.
