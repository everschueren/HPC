#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(stringr))
suppressMessages(library(sqldf))
suppressMessages(library(optparse))
suppressMessages(library(compiler))


Flatten.Fisher = function(x) {
  pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)  
}

Flatten.Resampling = function(){
  print("ERROR: NOT IMPLEMENTED YET")
  quit()
}

Flatten.getFlattenFun = function(fun_str){
  if(fun_str == "Resampling"){
    Flatten.Resampling    
  }else{
    Flatten.Fisher
  }
}

Flatten.main = function(data_file, output_file, method="Fisher"){
  data = read.delim(data_file, stringsAsFactors=F, check.names=F)
  flat_fun=Flatten.getFlattenFun(method)
  data_flat = sqldf("select uniprot_ac, count(*) as 'pep_count' from data group by uniprot_ac order by uniprot_ac asc")
  
  i = 4
  data_colnames = gsub("-|\\+","_",colnames(data))
  # print(data_colnames)
  colnames(data) = data_colnames
  
  while(i < ncol(data)){

    lfc_col = colnames(data)[i]
    p_col = colnames(data)[i+1]
    
    lfc_flat = sqldf(str_join("select uniprot_ac, avg(`",lfc_col,"`) as '",lfc_col,"_avg', stdev(`",lfc_col,"`) as '",lfc_col,"_std' from data where `",lfc_col,"` != 'NA' group by uniprot_ac order by uniprot_ac asc"))
    
    pvls = data[,c(2,i+1)]
    pvls_agg = aggregate(. ~ uniprot_ac, data=pvls, FUN=flat_fun)
    data_flat = merge(data_flat, lfc_flat, by=("uniprot_ac"))
    data_flat = merge(data_flat, pvls_agg, by=("uniprot_ac"))
    
    i=i+2
  }
  
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
Flatten.main <- cmpfun(Flatten.main)
Flatten.main(data_file=parsedArgs$data_file, output_file=parsedArgs$output_file, method=parsedArgs$method)  

# Flatten.main(data_file="~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-all-Ub/processed_repeats/Mock_v_WT_all_evidence_FLT_MAT_LIM.txt", output_file="~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-all-Ub/processed_repeats/Mock_v_WT_all_evidence_FLT_MAT_LIM_PRT.txt", method="Fisher")  
