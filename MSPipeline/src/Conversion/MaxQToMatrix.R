#! /usr/bin/Rscript --vanilla --default-packages=utils
## CONVERTS MAXQUANT OUTPUT TO MATRIX WITH INTENSITIES

library(sqldf)
library(reshape2)
library(stringr)
library(optparse)
library(compiler)

options(warn=-1)

max_ratios = function(vec){
  min_val = min(vec, na.rm=T)
  max_val = max(vec, na.rm=T)
  if(is.na(min_val) == F & is.infinite(min_val) == F){
    if(abs(log2(min_val)) >= abs(log2(max_val))){
      min_val
    }else{
      max_val
    }  
  }else{
    NA
  }
}


get_replicate_processing_fun = function(method="none"){
  if(method=="max_ratio"){
    fun = max_ratios
  }else if(method == "max"){
    fun = max
  }else if(method=="mean"){
    fun = mean
  }else if(method=="none"){
    fun = identity
  }
  fun
}

get_replicate_columns = function(keys, sample){
 cnames = as.vector(unlist(sqldf(str_join("select Code from keys where Sample = '",sample,"'")))) 
}

maxQ_to_matrix = function(data_file, index_file="", output_file="", replicate_processing="none", technical_replicates=2, maxq_value="Intensity",sequence_filter="none", uniqueness_filter=T){
  
  ## READ DATA
  data = read.delim(data_file, stringsAsFactors=F)
  print("SAMPLECODES IN DATAFILE:")
  print(unique(data$Raw.file))

  ## Extract all required columns and replace 'File' name by 'Code' 
  keys = read.delim(index_file, stringsAsFactors=F)
  
  ## CONVERT FROM LONG TO WIDE FORMAT
  fun = get_replicate_processing_fun(replicate_processing)
  tmp2 = dcast(tmp, uniprot_id+Sequence ~ Code,value.var=maxq_value,  fun.aggregate=fun, fill=1)
  data_keys = tmp2[,1:2]
  data_matrix = tmp2[,c(3:ncol(tmp2))]
  data_matrix[data_matrix==1]=1 
  ## FLATTEN OVER REPLICATES TO GET ONE VALUE ER BIOLOGICAL SAMPLE 
  data_matrix_flat = c()
  sample_keys = unique(keys$Sample)
  for(i in 1:length(sample_keys)){
    key = sample_keys[i]
    cnames = intersect(get_replicate_columns(keys,sample_keys[i]),colnames(data_matrix))
    data_matrix_sample = data_matrix[, cnames]
    col = apply(data_matrix_sample,1,fun)
    data_matrix_flat = cbind(data_matrix_flat, col)
    colnames(data_matrix_flat)[i] = key
  }
  
  data_matrix_flat = cbind(data_keys, data_matrix_flat)
  write.table(data_matrix_flat, file=output_file, eol="\n", quote=F, sep="\t", row.names=F, col.names=T)
}

## READ OPTIONS

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-i", "--index_file"),
              help="index file mapping raw files to sample codes"),
  make_option(c("-d", "--data_file"),
              help="data file containing maxquant output"),
  make_option(c("-o", "--output_file"),
              help="output file for converted matrix"),
  make_option(c("-r", "--replicate_processing"), default="none", 
              help="function to process replicates: none, max, mean"),
  make_option(c("-t", "--tech_reps"), default=2, type="integer", 
              help="number of technical replicates"),
  make_option(c("-m", "--maxq_value"), default="Intensity",
              help="String indicating whether to parse out the Intensity or Ratio_H_L column")
)


parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
maxQ_to_matrix <- cmpfun(maxQ_to_matrix)
maxQ_to_matrix(parsedArgs$data_file, parsedArgs$index_file, parsedArgs$output_file, parsedArgs$replicate_processing, parsedArgs$tech_reps, parsedArgs$maxq_value, parsedArgs$sequence_filter, parsedArgs$uniqueness_filter)  


