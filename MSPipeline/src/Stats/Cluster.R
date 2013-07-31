#! /usr/bin/Rscript --vanilla --default-packages=utils
INSTALL_DIR = Sys.getenv("MS_PIPELINE_PATH") 
source(paste(INSTALL_DIR,"src/Plotting/Heatmaps.R",sep=""))

suppressMessages(library(optparse))
suppressMessages(library(compiler))

options(warn=-1)

clusterIPs.main = function(data_file, output_file){
  
  results_MAT = read.delim(data_file, stringsAsFactors=F)
  results_baits = results_MAT[1,5:ncol(results_MAT)]
  
  results_MAT_clean = results_MAT[3:nrow(results_MAT),5:ncol(results_MAT)]
  colnames(results_MAT_clean) = paste(results_baits, colnames(results_baits))
  
  rownames(results_MAT_clean) = results_MAT[3:nrow(results_MAT),1]
  results_MAT_clean = data.matrix(results_MAT_clean)
  results_MAT_cor = cor(results_MAT_clean, use="pairwise.complete.obs", method="pearson")
  
  rowScale = 100 / (nrow(results_MAT_cor) * 2.2)
  pdf(output_file, width=7, height=7)
    heatmap.EV(results_MAT_cor, cexRow=rowScale)
  dev.off()
}

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--data_file"),
              help="data file containing matrix"),
  make_option(c("-o", "--output_file"),
              help="output file for cluster plot")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
clusterIPs.main(parsedArgs$data_file, parsedArgs$output_file)