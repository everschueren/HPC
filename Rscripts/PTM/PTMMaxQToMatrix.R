library(sqldf)
library(limma)
library(reshape2)
library(stringr)

scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
plotDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/output/"
processedDataDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/MAXQ_MATRICES/"
inputDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/input/"

##########################
## reformat to data matrix

MaxQToMatrix = function(maxQData){
  colnames(maxQData) = tolower(colnames(maxQData))
  evidence = sqldf("select raw_file, modified_sequence as 'mod_seq',Avg(intensity_l) as 'L',avg(intensity_h) as 'H' from maxQData group by raw_file, mod_seq")
  
  evidence.L = evidence[,c(1,2,3)]
  evidence.H = evidence[,c(1,2,4)]
  
  evidence.L.wide <- dcast(evidence.L, mod_seq  ~ raw_file, value.var="L")
  evidence.H.wide <- dcast(evidence.H, mod_seq  ~ raw_file, value.var="H")
  colnames(evidence.L.wide) = str_join(colnames(evidence.L.wide),"_L")
  colnames(evidence.H.wide)[2:5] = str_join(colnames(evidence.H.wide[2:5]),"_H")
  evidence.wide = cbind(evidence.H.wide, evidence.L.wide[,2:5])
  evidence.wide.NNA = evidence.wide[rowSums(is.na(evidence.wide))!=(ncol(evidence.wide)-1), ]
  evidence.wide.NNA  
}

maxQFiles = list.files(inputDir,pattern = "evidence.txt")
for(f in maxQFiles){
  maxQData <- read.delim(str_join(inputDir,f))
  evidence.wide = MaxQToMatrix(maxQData)
  write.table(evidence.wide,file=str_join(processedDataDir,f),sep="\t",eol="\n",quote=F,row.names=F)
}
