#!/bin/env/ Rscript

library(dplyr)
library(writexl)

parser <- argparse::ArgumentParser()
parser$add_argument("--gwas", type="character")
parser$add_argument("--marker", type="character")
parser$add_argument("--output", type="character")
args <- parser$parse_args()

if(!file.exists(args$gwas)) {
    stop("Error: GWAS input file does not exist!")
}

if(!file.exists(args$marker)) {
    stop("Error: marker input file does not exist!")
}

if(args$output == "") {
    stop("Error: No output filename has been given!")
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
    dplyr::mutate(anchored = "yes") %>%
    dplyr::filter(!is.na(FDR_BH))

lg_marker <- dplyr::inner_join(x = marker_position, y = lg_marker, by = "ID")

non_lg_marker <- data %>%
    dplyr::filter(!stringr::str_detect(CHR, "LG.*")) %>%
    dplyr::filter(!is.na(FDR_BH)) %>%
    dplyr::mutate(anchored = "no") %>%
    dplyr::mutate(Reference = CHR) %>%
    dplyr::mutate(Position = stringr::str_split(ID, "_") %>% purrr::map_chr(2)) %>%
    dplyr::mutate(across(c("UNADJ", "GC", "QQ", "BONF", "HOLM", "SIDAK_SS", "SIDAK_SD", "FDR_BH", "FDR_BY"), ~ formatC(.x, format = "e", digits = 3))) %>%
    dplyr::relocate(Reference, Position, ID, CHR)

output <- rbind(lg_marker, non_lg_marker) %>%
    dplyr::arrange(FDR_BH)

writexl::write_xlsx(output, args$output)
