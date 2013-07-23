#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(stringr))
suppressMessages(library(sqldf))
suppressMessages(library(optparse))
suppressMessages(library(compiler))


FlattenP.Fisher = function(x) {
  pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)  
}


Flatten = function(data_file, output_file, method="Fisher"){
  data = read.delim(data_file, stringsAsFactors=F, check.names=F)
  lfc_col = colnames(data)[4]
  data_flat = sqldf(str_join("select uniprot_ac, avg(`",lfc_col,"`) as 'LFC_avg', stdev(`",lfc_col,"`) as 'LFC_std', count(*) as 'pep_count' from data where `",lfc_col,"` != 'NA' group by uniprot_ac order by uniprot_ac asc"))
  if(method=="Fisher"){
    pvls = data[,c(2,5)] ## HARDCODED
    pvls_agg = aggregate(. ~ uniprot_ac, data=pvls, FUN=FlattenP.Fisher)
  }
  data_flat = cbind(data_flat, pvls_agg[,2])
  colnames(data_flat)[5] = 'P' 
  write.table(data_flat, file=output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
}

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--data_file"),
              help="Data file with limma-like output (Protein, Peptide, LFC, p-value)"),
  make_option(c("-o", "--output_file"),
              help="Output file"),
  make_option(c("-m", "--method"),
              help="Method to flatten peptides to proteins, options are: Fisher, Resampling")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
Flatten <- cmpfun(Flatten)
Flatten(data_file=parsedArgs$data_file, output_file=parsedArgs$output_file, method=parsedArgs$method)  

# Flatten(data_file="~/Projects/HPCKrogan/Scripts/MSPipeline/tests/apms_maxq/processed/vpr_timecourse_FLT_MAT_LIM.txt", output_file="~/Projects/HPCKrogan/Scripts/MSPipeline/tests/apms_maxq/processed/vpr_timecourse_FLT_MAT_LIM_AVG.txt", method="Fisher")  
