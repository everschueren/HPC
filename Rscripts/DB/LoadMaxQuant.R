library(RMySQL)
library(stringr)
library(sqldf)

## make db connection and get all data
## find a solution for the ssh tunnel that can't be run in R 
#ssh -L 3307:localhost:3306 everschueren@bluemoon.ucsf.edu

scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
processedDataDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/"

setwd(scriptDir)
source("DB/LoadUniprot.R")

connect = function(){
  con = dbConnect(MySQL(),host="127.0.0.1",dbname="HPCKrogan",user="everschueren",pass="hpckr0gan",port=3307)  
  con
}

disconnect = function(con){
  dbDisconnect(con)
}

LoadMaxQuantData = function(attributeList="*",selectCondition=""){
  con = connect()
  query = str_join("select ", attributeList ," from MS_Result_maxquant R join MS_Experiment E on E.ms_experiment_id = R.ms_experiment_id ", selectCondition)
  result = dbGetQuery(con,query)
  disconnect(con)
  result
}

LoadMaxQuantPeptides = function(){
  con = connect()
  query = "select peptide, proteins from MS_Result_maxquant group by peptide, proteins"
  result = dbGetQuery(con,query)
  disconnect(con)
  result
}

# uniprotData = write.table(uniprotData, file=str_join(processedDataDir,"uniprot_protein_descriptions.txt"), eol="\n", quote=F, sep="\t", row.names=F, col.names=T)
# uniprotData = read.delim(str_join(processedDataDir,"uniprot_protein_descriptions.txt"))

LoadMaxQuantPeptidesWithDescriptions = function(){
  
  maxQPeptides = LoadMaxQuantPeptides()
  uniprotData = LoadUniprotData()
  maxQDeconvoluted = c()
  
  ## loop over data frame, split and add 
  for(i in 1:nrow(maxQPeptides)){
    ## get key column
    d = maxQPeptides[i,]
    dkey = d[,"proteins"]
    
    ## deconvolute multiple matches in key column
    splitkeys = str_split(dkey, ';')[[1]]
    for(splitkey in splitkeys){
      maxQDeconvoluted = rbind(maxQDeconvoluted, c(d$peptide, splitkey))  
    }
  }
  
  maxQDeconvoluted = as.data.frame(maxQDeconvoluted)
  colnames(maxQDeconvoluted) = c("peptide","key")
  
  ## join to uniprot
  ## set up dummy for entries not found
  dummy = as.data.frame(matrix("",nrow=1,ncol=ncol(uniprotData)))
  colnames(dummy) = colnames(uniprotData)
  
  peptide_protein_table = sqldf("select * from maxQDeconvoluted M left join uniprotData U on M.key = U.uniprot_id group by M.peptide, M.key" , stringsAsFactors=F)
  uniprotData = write.table(peptide_protein_table, file=str_join(processedDataDir,"peptide_protein_descriptions.txt"), eol="\n", quote=F, sep="\t", row.names=F, col.names=T)
}


