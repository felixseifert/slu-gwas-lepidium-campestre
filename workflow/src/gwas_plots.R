#!/bin/env/ Rscript

# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Generation of Manhattan and QQ-plots for GWAS analysis

library(dplyr)

parser <- argparse::ArgumentParser()
parser$add_argument("--gwas", type="character")
parser$add_argument("--marker", type="character")
parser$add_argument("--output_prefix", type="character")
args <- parser$parse_args()

if(!file.exists(args$gwas)) {
    stop("Error: GWAS input file does not exist!")
}

if(!file.exists(args$marker)) {
    stop("Error: marker input file does not exist!")
}

if(args$output_prefix == "") {
    stop("Error: No output filename prefix has been given!")
}

data <- read.table(
    args$gwas,
    header = TRUE
)
colnames(data) <- c("CHR", "ID", "UNADJ", "GC", "QQ", "BONF", "HOLM", "SIDAK_SS", "SIDAK_SD", "FDR_BH", "FDR_BY")

marker_position <- read.csv(args$marker, sep="\t") %>%
    dplyr::select(Reference, Position, ID)

lg_marker <- data %>%
    dplyr::filter(stringr::str_detect(CHR, "LG.*")) %>%
    dplyr::mutate(CHR = stringr::str_replace(CHR, 'LG', '')) %>%
    dplyr::filter(!is.na(FDR_BH))

lg_marker <- dplyr::inner_join(x = marker_position, y = lg_marker, by = "ID")

lg_marker$CHR <- as.numeric(factor(lg_marker$CHR))

png(paste0(args$output_prefix, "_manhattan.png"))
qqman::manhattan(lg_marker, chr="CHR", bp="Position", p="FDR_BH", snp="ID")
dev.off()

png(paste0(args$output_prefix, "_qq.png"))
qqman::qq(lg_marker$FDR_BH)
dev.off()
