#!/usr/bin/env Rscript

# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# Generates pophelper plots

library("optparse")
library("stringr")

option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="clumpp formatted input file(s) for pophelper, multiple files separated by comma", metavar="character"),
	make_option(c("-o", "--outdir"), type="character", default=NULL, help="output directory for plot", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$input) || is.null(opt$outdir)) {
	print_help(opt_parser)
	stop("The input and output filename need to be given as parameters.", call.=FALSE)
}

library("pophelper")

input_files = str_split(opt$input, ",")[[1]]
data = readQ(toString(input_files[1]))
for(index in 2:length(input_files)) {
    data = joinQ(data, readQ(input_files[index]))        
}

cluster_number = function(x) attr(x, "k")
spnames = paste0("K=",sapply(data, cluster_number))

output_directory = paste(opt$outdir, "/", sep="")

plotQ(data, basesize=5, dpi=300, exportpath=output_directory, imgoutput="join", imgtype="png", returnplot=T, sharedindlab=F, sortind="all", splab=spnames, showlegend=T, legendtextsize=3, legendmargin=c(1,1,1,0), showyaxis=T, showticks=T, height=2, width=10)
plotQ(data, basesize=5, dpi=300, exportpath=output_directory, imgoutput="join", imgtype="pdf", returnplot=T, sharedindlab=F, sortind="all", splab=spnames, showlegend=T, legendtextsize=3, legendmargin=c(1,1,1,0), showyaxis=T, showticks=T, height=2, width=10)

plotQMultiline(data, basesize=5, dpi=300, exportpath=output_directory, imgtype="png", indlabsize=8, legendmargin=c(1,1,1,0), legendtextsize=10, returnplot=T, showlegend=T, sortind="all")
plotQMultiline(data, basesize=5, dpi=300, exportpath=output_directory, imgtype="pdf", indlabsize=8, legendmargin=c(1,1,1,0), legendtextsize=10, returnplot=T, showlegend=T, sortind="all")
