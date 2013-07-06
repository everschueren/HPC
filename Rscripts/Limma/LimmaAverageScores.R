#! /usr/bin/Rscript --vanilla --default-packages=utils

library(stringr)
library(sqldf)
library(optparse)
library(compiler)

select_and_flatten = function(limma_result_file, p_value, num_peptides=1, output_dir){
  limma_result = read.delim(limma_result_file)
  desc = limma_result[,1:5]
  data = limma_result[,6:ncol(limma_result)]
  bname = basename(limma_result_file)
  bname = substr(bname,1,nchar(bname)-4)
  ## data is in format (logFC+p-val)1 (logFC+p-val)2 (logFC+p-val)3 ...
  for(i in 1:(ncol(data)/2)){
    start = (i*2)-1
    end = i*2
    tmp  = cbind(desc,data[,start:end])
    cname=colnames(data)[start]
    name=substr(cname,1,nchar(cname)-6)
    colnames(tmp)[6:7] = c("logFC","PValue")
    contrast = sqldf(str_join("select uniprot_ac, description, gene_name, uniprot_id, avg(logFC) as logFC_avg, stdev(logFC) as logFC_std ,count(*) as 'num_peptides' from tmp where logFC>0 and PValue <= ", p_value, " group by uniprot_id having count(*) >=  ", num_peptides, " order by count(*) desc"))
    write.table(contrast, file=str_join(output_dir,"/",bname,"_",name,"_selected_flat.txt"), row.names=F, col.names=T,quote=F,sep="\t",eol="\n")
  }
}


option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-l", "--limma_file"),
              help="result file containing output values from Limma"),
  make_option(c("-p", "--p_value"), type="double", default=0.05,
              help="p-value to select below on the peptide level before flattening"),
  make_option(c("-n", "--num_peptides"), type="integer", default=1,
              help="number of peptides a protein should have for being selected"),
  make_option(c("-d", "--output_dir"),
              help="output dir for flattened contrasts")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))

## MAIN
select_and_flatten_C <- cmpfun(select_and_flatten)
select_and_flatten_C(parsedArgs$limma_file, parsedArgs$p_value, parsedArgs$num_peptides, parsedArgs$output_dir)

 
