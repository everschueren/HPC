#! /usr/bin/Rscript --vanilla --default-packages=utils

INSTALL_DIR = Sys.getenv("MS_PIPELINE_PATH") 
#source(paste(INSTALL_DIR,"src/Plotting/Heatmaps.R",sep=""))

suppressMessages(library(optparse))
suppressMessages(library(compiler))
suppressMessages(library(stats))
suppressMessages(library(grDevices))
suppressMessages(library(graphics))
suppressMessages(library(pheatmap))

options(warn=-1)

clusterIPs.main = function(data_file, output_file, font_scale){
  
  results_MAT = read.delim(data_file, stringsAsFactors=F)
  ## strip out unnamed Ips 
  results_MAT = results_MAT[,results_MAT[2,]!="" & is.na(results_MAT[2,])==F]
  results_baits = results_MAT[1,5:ncol(results_MAT)]
  
  results_MAT_clean = results_MAT[3:nrow(results_MAT),5:ncol(results_MAT)]
  colnames(results_MAT_clean) = paste(results_baits, colnames(results_baits))
  
  rownames(results_MAT_clean) = results_MAT[3:nrow(results_MAT),1]

  results_MAT_clean = data.matrix(results_MAT_clean)
  results_MAT_cor = cor(results_MAT_clean, use="pairwise.complete.obs", method="pearson")
  
  #rowScale = font_scale / (nrow(results_MAT_cor))
  #rowScale = (.9406 -.00508*nrow(results_MAT_cor))	#basic linear modelsolution from other sizes
  
  pheatmap(results_MAT_cor, cluster_rows=T, cluster_cols=T, scale="none",fontsize_row=font_scale,fontsize_col=font_scale, cellwidth=font_scale, cellheight=font_scale, border_color=NA, filename=output_file)

}

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--data_file"),
              help="data file containing matrix"),
  make_option(c("-o", "--output_file"),
              help="output file for cluster plot"),
  make_option(c("-s", "--font_scale"), default=40,
              help="scaling factor for fonts on the rows and columns")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))

TEST=F
if(TEST){
  parsedArgs$data_file = "~/Projects/HPCKrogan/Data/HHV8/data/iSLK/data/processed/iSLK_data_wKEYS_NoC_MAT.txt"
  data_file=parsedArgs$data_file
  parsedArgs$output_file = "~/Projects/HPCKrogan/Data/HHV8/data/iSLK/data/processed/iSLK_data_wKEYS_NoC_MAT_CLUSTERED.pdf"
  output_file = parsedArgs$output_file
}

clusterIPs.main(parsedArgs$data_file, parsedArgs$output_file, parsedArgs$font_scale)

