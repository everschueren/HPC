#! /usr/bin/Rscript --vanilla --default-packages=utils
## CONVERTS MAXQUANT OUTPUT TO MATRIX WITH INTENSITIES

suppressMessages(library(sqldf))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(compiler))
suppressMessages(library(limma))
suppressMessages(library(stats))

options(warn=-1)

MaxQToWide = function(data=tmp, id_idx=c(1,2), sample_idx=3, value_idx=4, fun=mean){
  print(">> BUILDING MATRIX")
  ## get unique id combinations and create a unique key
  id_col_names = colnames(data)[id_idx]
  sample_col_name = colnames(data)[sample_idx]
  unique_id_combinations = unique(data[,id_idx])
  unique_id_combinations = cbind(1:nrow(unique_id_combinations),unique_id_combinations)
  colnames(unique_id_combinations)[1] = "unique_id"
  
  ## merge keys with the data
  data_long = merge(unique_id_combinations, data, by=id_col_names, all.x=T)[,(length(id_idx)+1):(length(id_idx)+3)]
  sample_names = unique(data_long[,sample_col_name])
  
  wide = NULL
  
  for(i in 1:length(sample_names)){
    sample_name = sample_names[i]
    print(sprintf("adding sample %s",sample_name))
    sub_set = data_long[data_long[,sample_col_name]==sample_name,]
    tmp = unique_id_combinations[,"unique_id", drop=F]
    data_long_sample = merge(tmp,sub_set, by="unique_id", all.x=T)
    data_long_sample = aggregate(data_long_sample$Ratio_H_L, by=list(unique_id=data_long_sample$unique_id), FUN=fun)
    
    if(is.null(wide)){
      wide = data_long_sample
    }else{
      wide = cbind(wide, data_long_sample$x)
    }
    colnames(wide)[ncol(wide)] = sample_name
  }
  wide = merge(unique_id_combinations, wide, by=c("unique_id"))
  wide = wide[,2:ncol(wide)]
  wide
}
MaxQToWide = cmpfun(MaxQToWide)

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
    fun = function(x) max(x, na.rm=T)
  }else if(method=="mean"){
    fun = function(x) mean(x, na.rm=T)
  }else if(method=="none"){
    fun = identity
  }
  fun
}

get_replicate_columns = function(keys, sample){
 cnames = as.vector(unlist(sqldf(str_join("select Code from keys where Sample = '",sample,"'")))) 
}

maxQ_to_matrix = function(data_file, index_file="", output_file="", replicate_processing="none", maxq_value="Ratio_H_L", normalization="none"){
  
  ## READ DATA
  data = read.delim(data_file, stringsAsFactors=F)
  print("SAMPLECODES IN DATAFILE:")
  print(unique(data$Raw.file))
  
  ## Extract all required columns and replace 'File' name by 'Code' 
  keys = read.delim(index_file, stringsAsFactors=F)
  #tmp = sqldf(str_join("select  Proteins as 'uniprot_id', Modified_Sequence as 'Sequence', Code, ",maxq_value," from data D join keys K on K.File = D.Raw_File where ",maxq_value," != 'NA'"))
  tmp = sqldf(str_join("select  Proteins as 'uniprot_ac', Modified_sequence as 'Sequence', Code, ",maxq_value," from data D join keys K on K.File = D.Raw_file"))
  print(colnames(tmp))
  print(nrow(tmp))
  
  ## CONVERT FROM LONG TO WIDE FORMAT
  fun = get_replicate_processing_fun(replicate_processing)
  # tmp2 = dcast(tmp, uniprot_id+Sequence ~ Code,value.var=maxq_value,  fun.aggregate=fun, fill=1)
  
  ## dcast crashes on big tables, therefore we use our won function
  #tmp2 = dcast(tmp, uniprot_ac+Sequence ~ Code,value.var=maxq_value,  fun.aggregate=fun)
  tmp2 = MaxQToWide(tmp, id_idx=c(1,2), sample_idx=3, value_idx=4, fun=fun)
  
  data_keys = tmp2[,1:2]
  data_matrix = tmp2[,c(3:ncol(tmp2))]
  #data_matrix[data_matrix==1]=1 

  ## FLATTEN OVER TECH REPLICATES TO GET ONE VALUE PER BIOLOGICAL SAMPLE 
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
  
  ## NORMALIZE
  data_matrix_flat[is.infinite(data_matrix_flat) | is.nan(data_matrix_flat)] = NA
  if(normalization != "none"){
    data_matrix_flat = normalizeBetweenArrays(data_matrix_flat, method=normalization)
  }
  ## CONVERT TO LOGS
  data_matrix_flat = log2(data_matrix_flat)
  output = cbind(data_keys, data_matrix_flat)
  write.table(output, file=output_file, eol="\n", quote=F, sep="\t", row.names=F, col.names=T)
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
  make_option(c("-m", "--maxq_value"), default="Intensity",
              help="String indicating whether to use Intensities or Ratio_H_L column to compute ratios"),
  make_option(c("-n", "--normalization"), default="none",
              help="String indicating which normalization method to use")
)


parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
maxQ_to_matrix <- cmpfun(maxQ_to_matrix)
maxQ_to_matrix(parsedArgs$data_file, parsedArgs$index_file, parsedArgs$output_file, parsedArgs$replicate_processing, parsedArgs$maxq_value, parsedArgs$normalization)  

# maxQ_to_matrix(index_file='~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-Ub/input/Mock_v_WT_keys.txt', data_file='~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-Ub/input/Mock_v_WT_evidence.txt', output_file='~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-Ub/processed_repeats/Mock_v_WT_evidence_FLT_MAT.txt', replicate_processing='mean', maxq_value='Ratio_H_L', normalization='scale')
