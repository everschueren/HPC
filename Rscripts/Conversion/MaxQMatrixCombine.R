#! /usr/bin/Rscript --vanilla --default-packages=utils
## CONVERTS MAXQUANT OUTPUT TO MATRIX WITH INTENSITIES

library(sqldf)
library(reshape2)
library(stringr)
library(optparse)
library(compiler)

maxQ_matrix_combine = function(index_file, output_file, join="inner"){
  index = read.delim(index_file, stringsAsFactors=F)
  tmp = NULL
  for(f in index$Files){
    file = read.delim(f, stringsAsFactors=F)
    f_col_names = colnames(file)[3:ncol(file)]
    f_col_names_select = str_join("F.",f_col_names,collapse=",")
    if(is.null(tmp)){
      tmp = file
    }else{
      #tmp = sqldf(str_join("select T.*, ",f_col_names_select," from tmp T ",join," join file F on (T.Sequence = F.Sequence and T.uniprot_id = F.uniprot_id) "))
      #print(colnames(tmp))
      #print(colnames(file))
      
      tmp = merge(tmp, file, by=c("Sequence", "uniprot_id"), all.x=T, all.y=T)
    }
  }
  tmp[is.na(tmp)]=1
  write.table(tmp, file=output_file, eol="\n", quote=F, sep="\t", row.names=F, col.names=T) 
}

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-i", "--index_file"),
              help="index file listing processex MaxQ to Matrix files"),
  make_option(c("-o", "--output_file"),
              help="output file for Joined Matrix Files on Sequence and Uniprot_id"),
  make_option(c("-j", "--join"),default = "inner",
              help="join type to combine matrix files: inner/left")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
maxQ_matrix_combine_C <- cmpfun(maxQ_matrix_combine)
maxQ_matrix_combine_C(parsedArgs$index_file, parsedArgs$output_file, parsedArgs$join)