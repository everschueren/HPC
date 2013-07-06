#! /usr/bin/Rscript --vanilla --default-packages=utils
## CONVERTS MAXQUANT OUTPUT TO MATRIX WITH INTENSITIES

library(sqldf)
library(reshape2)
library(stringr)
library(optparse)
library(compiler)

options(warn=-1)

explode_pairs = function(pair){
  value = pair[[1]]
  splitkeys = str_split(pair[[2]], ';')[[1]]
  newpair = c()
  for(splitkey in splitkeys){
    newpair = c(newpair, c(value, splitkey))  
  }
  newpair
}

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

filter_peptides_by_modification = function(data, filter){
  print(str_join("FILTER: ",filter))
  if(filter == "UBI"){
    pattern = "'%(gl)%'"
  }else if(filter == "ACE"){
    pattern = "'%(ac)%'"
  }else{
    pattern = "'%'"
  }
  sqldf(str_join("select * from data where Sequence like ",pattern))
}

filter_peptides_by_contaminants = function(data){
  print(str_join("FILTER: ","CON__,REV__"))
  sqldf(str_join("select * from data where uniprot_id not like '%CON__%' and uniprot_id not like '%REV__%'"))
}

filter_peptides_by_uniqueness = function(data){
  sqldf(str_join("select * from data where uniprot_id not like '%;%'"))
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
  tmp = sqldf(str_join("select  Proteins as 'uniprot_id', Modified_Sequence as 'Sequence', Code, ",maxq_value," from data D join keys K on K.File = D.Raw_File where ",maxq_value," != 'NA'"))
  ## FILTER ALL UNWANTED SEQUENCES AND IDENTIFIED PROTEINS
  tmp = filter_peptides_by_contaminants(tmp)
  tmp = filter_peptides_by_modification(tmp, sequence_filter)
  if(uniqueness_filter){
    print(str_join("FILTER: ","UNIQUE_PEPTIDES"))
    tmp = filter_peptides_by_uniqueness(tmp)  
  }
  
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


## PART TO EXPLODE PAIRS: 
#   maxq_unique_assigned = c()
#   
#   ## loop over data frame, split and add 
#   for(p in 1:nrow(keys)){
#       maxq_unique_assigned = c(maxq_unique_assigned, explode_pairs(keys[p,]))  
#   }
#   keys = as.data.frame(matrix(maxq_unique_assigned, ncol=2, nrow=length(maxq_unique_assigned)/2, byrow=T))
#   colnames(keys) = c("Sequence", "uniprot_id")
# result = sqldf("select K.Sequence, K.uniprot_id, D.* from data_keys K join data_matrix_flat D on K.Sequence = D.Sequence where K.uniprot_id not like '%CON__%' and K.uniprot_id not like '%REV__%'")


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
              help="String indicating whether to parse out the Intensity or Ratio_H_L column"),
  make_option(c("-f", "--sequence_filter"), default="none",
              help="String indicating a pattern to select for specific modified sequences like (gl)=>ubiquilyated"),
  make_option(c("-u", "--uniqueness_filter"), default=TRUE,
              help="boolean indicating a whether to only select uniquely assigned peptides")
)


## TEST ##############################################
TEST=FALSE

if(TEST){
  TESTDIR="~/Projects/HPCKrogan/Data/Mtb-proteomics/TJP1-4/"
  PREFIX="082312-TJP-1-4-Ub"
  NORMALIZATION="cyclicloess"
  REPLICACOMBINATION="max_ratio"
  REPLICAS=2
  PVALUE=0.05
  MAXQCOLUMN="Ratio_H_L"
  FILTER="UBI"
  test_args = c("-i",str_join("input/",PREFIX,"-files.txt"), "-d", str_join("input/",PREFIX,"-evidence.txt"), "-o",str_join("processed/",PREFIX,"-matrix-max.txt"), "-r", REPLICACOMBINATION, "-t", REPLICAS, "-m", MAXQCOLUMN, "-f",FILTER)
  setwd(TESTDIR)
  parsedArgs = parse_args(OptionParser(option_list = option_list), args = test_args)  
  debug(maxQ_to_matrix)
  maxQ_to_matrix(parsedArgs$data_file, parsedArgs$index_file, parsedArgs$output_file, parsedArgs$replicate_processing, parsedArgs$tech_reps, parsedArgs$maxq_value, parsedArgs$sequence_filter)
}else{
  parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
  maxQ_to_matrix_C <- cmpfun(maxQ_to_matrix)
  maxQ_to_matrix_C(parsedArgs$data_file, parsedArgs$index_file, parsedArgs$output_file, parsedArgs$replicate_processing, parsedArgs$tech_reps, parsedArgs$maxq_value, parsedArgs$sequence_filter, parsedArgs$uniqueness_filter)  
}


