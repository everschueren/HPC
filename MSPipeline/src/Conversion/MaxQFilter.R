#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(optparse))
suppressMessages(library(compiler))

options(warn=-1)

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


explode_pairs = function(pair){
  value = pair[[1]]
  splitkeys = str_split(pair[[2]], ';')[[1]]
  newpair = c()
  for(splitkey in splitkeys){
    newpair = c(newpair, c(value, splitkey))  
  }
  newpair
}

filter_peptides_by_modification = function(data, filter){
  print(paste("FILTER MODIFICATIONS :",filter))
  
  if(filter == "UBI"){
    pattern = '\\(gl\\)'
  }else if(filter == "ACE"){
    pattern = '\\(ac\\)'
  }else{
    pattern = '*'
  }
  keys = grep(pattern, data$Modified.sequence) 
  tmp = data[keys,]
  tmp
}

filter_peptides_by_contaminants = function(data){
  print("FILTER CONTAMINANTS")
  con_exp = 'CON__|REV__'  
  keys = grep(con_exp, data$Proteins, invert=T) 
  tmp = data[keys,]
  tmp
}


filter_peptides_by_uniqueness = function(data){
  print("FILTER UNIQUE PEPTIDES")
  unique_exp = ';'  
  keys = grep(unique_exp, data$Proteins, invert=T) 
  tmp = data[keys,]
  tmp
}

filter_peptides_by_oxidation = function(data){
  print("FILTER OXIDATIONS")
  ox_exp = '\\(ox\\)'  
  data$Modified.sequence = gsub(ox_exp, data$Modified.sequence, replacement='')
  data
}

filterData = function(data_file, output_file, contaminants_filter=T, uniqueness_filter=T, oxidations_filter=T, sequence_filter="none"){
  data = read.delim(data_file, stringsAsFactors=F)
  if(contaminants_filter){
    data = filter_peptides_by_contaminants(data)
  }
  if(uniqueness_filter){
    data = filter_peptides_by_uniqueness(data)
  }
  if(oxidations_filter){
    data = filter_peptides_by_oxidation(data)
  }
  data = filter_peptides_by_modification(data, sequence_filter)
  write.table(data, file=output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
} 

## READ OPTIONS

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--data_file"),
              help="data file containing maxquant output"),
  make_option(c("-o", "--output_file"),
              help="output file for filtered input"),
  make_option(c("-s", "--sequence_filter"), default="none",
              help="String indicating a pattern to select for specific modified sequences like (gl)=>ubiquilyated"),
  make_option(c("-u", "--uniqueness_filter"), default=TRUE,
              help="boolean indicating a whether to only select uniquely assigned peptides"),
  make_option(c("-c", "--contaminants_filter"), default=TRUE,
              help="boolean indicating a whether to filter out REV__/CON__"),
  make_option(c("-x", "--oxidations_filter"), default=TRUE,
              help="boolean indicating a whether to strip all (ox)")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
filterData <- cmpfun(filterData)
filter_peptides_by_modification <- cmpfun(filter_peptides_by_modification)
filter_peptides_by_oxidation <- cmpfun(filter_peptides_by_oxidation)
filter_peptides_by_uniqueness <- cmpfun(filter_peptides_by_uniqueness)
filter_peptides_by_contaminants <- cmpfun(filter_peptides_by_contaminants)
filterData(parsedArgs$data_file, parsedArgs$output_file, parsedArgs$contaminants_filter, parsedArgs$uniqueness_filter,  parsedArgs$oxidations_filter, parsedArgs$sequence_filter)  



