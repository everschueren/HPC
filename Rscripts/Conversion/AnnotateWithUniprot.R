library(sqldf)
library(stringr)
scriptdir = "~/Projects/HPCKrogan/Scripts/Rscripts/"

annotate_with_uniprot = function(data, key="uniprot_id", organism){
  organism_split = unlist(str_split(organism, "-"))
  #print(organism_split)
  Uniprot = NULL
  for(org in organism_split){
    tmp = read.delim2(str_join(scriptdir,"Conversion/Uniprot/uniprot_protein_descriptions_",org,".txt"), stringsAsFactors=F, quote="")    
    if(is.null(Uniprot)){
      Uniprot = as.data.frame(tmp)
    }else{
      Uniprot = rbind(Uniprot, tmp)  
    }
  }
  #print(nrow(Uniprot))
  sqldf(str_join("select D.",key," , U.Entry as 'uniprot_ac', U.Entry_name as 'uniprot_id', U.Gene_names as 'gene_name', U.Protein_names as 'description', U.Length as 'length', D.* from data D left join  Uniprot U on U.Entry=D.",key))
}