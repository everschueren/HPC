#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(sqldf))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(compiler))

annotate_with_uniprot = function(data, species="HUMAN", key="uniprot_ac", output_file=NULL, uniprot_dir="~/Projects/HPCKrogan/Scripts/MSPipeline/files/"){
  species_split = unlist(str_split(species, "-"))
  #print(organism_split)
  Uniprot = NULL
  for(org in species_split){
    tmp = read.delim2(str_join(uniprot_dir,"uniprot_protein_descriptions_",org,".txt"), stringsAsFactors=F, quote="")    
    if(is.null(Uniprot)){
      Uniprot = as.data.frame(tmp)
    }else{
      Uniprot = rbind(Uniprot, tmp)  
    }
  }
  #print(nrow(Uniprot))
  result = sqldf(str_join("select D.* , U.Entry as 'uniprot_ac', U.Entry_name as 'uniprot_id', U.Gene_names as 'gene_name', U.Protein_names as 'description', U.Length as 'length' from data D left join  Uniprot U on U.Entry=D.",key))
  if(is.null(output_file)){
    result
  }else{
    write.table(result, file=output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  }
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
  make_option(c("-u", "--uniprot_dir"),
              help="Organisms for which descriptions are available in the uniprot_dir"),
  make_option(c("-s", "--species"),
              help="Species for which descriptions are available in the uniprot_dir"),
  make_option(c("-k", "--key"),
              help="column containing the key in the data file")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
annotate_with_uniprot <- cmpfun(annotate_with_uniprot)
data = read.delim(parsedArgs$data_file, stringsAsFactors=F)
annotate_with_uniprot(data, species=parsedArgs$species, key=parsedArgs$key, output_file=parsedArgs$output_file, uniprot_dir=parsedArgs$uniprot_dir)  
