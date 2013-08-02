suppressMessages(library(sqldf))
suppressMessages(library(stringr))


annotate_with_uniprot = function(data, species="HUMAN", key="uniprot_ac", output_file=NULL, uniprot_dir="~/Projects/HPCKrogan/Scripts/MSPipeline/files/"){
  
  # print(colnames(data))
  data = as.data.frame(data)

  species_split = unlist(str_split(species, "-"))
  
  #print(organism_split)
  Uniprot = NULL
  for(org in species_split){
    print(paste("LOADING",org),sep="\t")
    tmp = read.delim2(str_join(uniprot_dir,"uniprot_protein_descriptions_",org,".txt"), stringsAsFactors=F, quote="")    
    if(is.null(Uniprot)){
      Uniprot = as.data.frame(tmp)
    }else{
      Uniprot = rbind(Uniprot, tmp)  
    }
  }
  # print(nrow(Uniprot))
  # print(nrow(data))
  result =  sqldf(str_join("select D.*, U.Entry as 'uniprot_ac', U.Entry_name as 'uniprot_id', U.Gene_names as 'gene_name', U.Protein_names as 'description', U.Length as 'length', substr(U.Entry_name,charindex('_',U.Entry_name)+1,length(U.Entry_name)) as 'species', substr(U.Entry_name,0,charindex('_',U.Entry_name,0)) as 'display_name' from data D left join  Uniprot U on U.Entry=D.",key))
  if(is.null(output_file)){
    result
  }else{
    write.table(result, file=output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  }
}