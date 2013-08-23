#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(reshape2))
suppressMessages(library(optparse))
suppressMessages(library(compiler))

options(warn=-1)

Report.cleanMatrix = function(data){
  ## get rid of unassigned columns
  bait_row = data[1,]
  data_tmp = data[,bait_row != ""]
  ## get rid of specificity exclusion row and meta-columns
  data_tmp = data_tmp[c(-1:-2),c(-2:-4)]
  colnames(data_tmp)[1] = "Preys"
  data_tmp
}

Report.IPoccurences = function(data){
  
  res = c()
  r_names = colnames(data[2:ncol(data)])
  for(i in 1:nrow(data)){
    r_key = data[i,1]
    r_data = data[i,2:ncol(data)]
    r_names_nonzero = r_names[r_data!=0]
    r_names_nonzero_coll = as.vector(paste(r_names_nonzero, collapse=","))
    res = rbind(res, as.vector(c(r_key, r_names_nonzero_coll)))  
  } 
  colnames(res) = c("PREY", "IPs")
  res
}
Report.IPoccurences = cmpfun(Report.IPoccurences)

Report.main = function(data_file, output_file){
  data = read.delim(data_file,  stringsAsFactors=F, header=T)
  ## convert into intermediate format
  data_tmp = Report.cleanMatrix(data)
  IP_occurences = Report.IPoccurences(data_tmp)
  write.table(IP_occurences, file=output_file, quote=F, row.names=F, sep="\t")
}


option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--data_file"),
              help="data file containing maxquant output"),
  make_option(c("-o", "--output_file"),
              help="output file for converted matrix")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
Report.main(parsedArgs$data_file, parsedArgs$output_file)

# Report.main(data_file="~/Projects/HPCKrogan/Data/HCV/Data/processed_293_andy/HCV-293T-Andy-results_wKEYS_NoC_MAT.txt", output_file="~/Projects/HPCKrogan/Scripts/MSPipeline/tests/report/HCV-293T-Andy-results_wMT_Report.txt")