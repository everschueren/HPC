library(sqldf)
library(optparse)
library(compiler)
library(limma)

get_timepoints_significant_proteins = function(significant_set_file, data_set_file, output_file, convert_to_logs=T, normalize_method="quantile"){
  significant_set =  read.delim(significant_set_file)
  data_set = read.delim(data_set_file)
  significant_data_set = sqldf("select S.uniprot_ac, D.* from significant_set S join data_set D on D.uniprot_id = S.uniprot_id ")
  if(convert_to_logs){
    tmp = log2(significant_data_set[,4:ncol(significant_data_set)]+1)
  }else{
    tmp = significant_data_set[,4:ncol(significant_data_set)]
  }
  if(normalize_method != "none"){
    tmp_colnames = colnames(tmp)
    tmp  = as.data.frame(normalizeBetweenArrays(as.matrix(tmp), method=normalize_method))
    colnames(tmp) = tmp_colnames
  }
  significant_data_set = cbind(significant_data_set[,1],tmp)
  colnames(significant_data_set)[1] = "uniprot_ac"
  significant_data_set_flat = aggregate(. ~ uniprot_ac,data=significant_data_set, FUN=mean)
  significant_data_set_flat_infected = significant_data_set_flat[,c(1,2,5,8,11)]
  plot_set = significant_data_set_flat_infected
  write.table(plot_set, file=output_file, quote=F,row.names=F,col.names=T, eol="\n", sep="\t")
}

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-s", "--significant_set_file"),
              help="result file containing averaged values after significant set computed by Limma"),
  make_option(c("-d", "--data_set_file"), 
              help="data file containing the reformated maxquant output with all timepoints"),
  make_option(c("-o", "--output_file"),
              help="output file for all timepoints"),
  make_option(c("-l", "--convert_to_logs"), default=T,
              help="convert maxquant values to log2-values"),
  make_option(c("-n", "--normalize_method"), default="none",
              help="Normalization method: quantile, cyclicloess, none")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))

## MAIN
get_timepoints_significant_proteins_C <- cmpfun(get_timepoints_significant_proteins)
get_timepoints_significant_proteins_C(parsedArgs$significant_set_file, parsedArgs$data_set_file, parsedArgs$output_file, parsedArgs$convert_to_logs, normalize_method=parsedArgs$normalize_method)
