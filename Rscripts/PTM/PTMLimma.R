library(stringr)
library(limma)
library(sqldf)

scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
plotDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/output/"
processedDataDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/LIMMA_NORMALIZED/"

setwd(scriptDir)
source("DB/LoadMaxQuant.R")
source("Plotting/Aux.R")
source("PTM/PTMAux.R")

# target = "VIF"
# target = "VPR"
target = "VPU"

#################################
## get all data in format L/H H/L

# IF NOT CACHE THEN LOAD FROM DB AND UPDATE TXT FILES ON DISK
# IF CACHE THEN READ FROM DISK
# IF NORMALIZED THEN LOAD avg_log_ratio_normalized FROM DB ELSE LOAD avg_log_ratio

CACHE=T
NORMALIZED=T

if(CACHE == F){
  MaxQData1 = LoadMaxQuantData()
  MaxQData2 = MaxQData1
  
  if(NORMALIZED){
    normalized_name = "_NORM_"
    dataset = sqldf(str_join("select R1.ms_internal_id, R1.proteins, R1.peptide, avg(R1.avg_log_ratio_normalized) as 'GFP_H_TARGET_L', avg(R2.avg_log_ratio_normalized) as 'TARGET_H_HIV_L' from MaxQData1 R1 left join MaxQData2 R2 on R1.`peptide` = R2.`peptide` where R1.`silac_heavy` = 'GFP' and R1.`silac_light` = '",target,"' and R2.`silac_heavy` = '",target,"' and R2.`silac_light` = 'GFP' and R1.ms_internal_id = R2.ms_internal_id and R1.peptide like '%(gl)%' group by ms_internal_id, peptide"))
  }else{
    dataset = sqldf(str_join("select R1.ms_internal_id, R1.proteins, R1.peptide, avg(R1.avg_log_ratio) as 'GFP_H_TARGET_L', avg(R2.avg_log_ratio) as 'TARGET_H_HIV_L' from MaxQData1 R1 left join MaxQData2 R2 on R1.`peptide` = R2.`peptide` where R1.`silac_heavy` = 'GFP' and R1.`silac_light` = '",target,"' and R2.`silac_heavy` = '",target,"' and R2.`silac_light` = 'GFP' and R1.ms_internal_id = R2.ms_internal_id and R1.peptide like '%(gl)%' group by ms_internal_id, peptide"))
    normalized_name = ""
  }
  setwd(processedDataDir)
  write.table(dataset, file=str_join(processedDataDir,target,"_UBI",normalized_name,".txt"), eol="\n", sep="\t", row.names=F, col.names=T, quote=F)
}else if(CACHE){
  setwd(processedDataDir)
  dataset = read.delim(str_join(processedDataDir,target,"_UBI",normalized_name,".txt"))
}

peptide_protein_descriptions <- read.delim(str_join(processedDataDir,"peptide_protein_descriptions.txt"))

#################################################################
## MERGE ALL DATA AS LEFT JOINS (NA for not corresponding fields)
## merge condition 1 and 2
dataset_merged = merge(dataset[dataset$ms_internal_id==str_join("PTM_UBI_",target,"_UNTREATED"),], dataset[dataset$ms_internal_id==str_join("PTM_UBI_",target,"_MG132"),], by="peptide", all.x=T, all.y=T)
dataset_merged = dataset_merged[, c(1,4,5,8,9)] 
colnames(dataset_merged)[2:5] = c(str_join("GFP_H_",target,"_L_UNTREATED"),str_join(target,"_H_GFP_L_UNTREATED"),str_join("GFP_H_",target,"_L_MG132"),str_join(target,"_H_GFP_L_MG132"))
## merge condition 3
dataset_merged = merge(dataset_merged, dataset[dataset$ms_internal_id==str_join("PTM_UBI_",target,"_IFN"),], by="peptide", all.x=T, all.y=T)
dataset_merged = dataset_merged[, c(1:5,8,9)] 
colnames(dataset_merged)[6:7] = c(str_join("GFP_H_",target,"_L_IFN"),str_join(target,"_H_GFP_L_IFN"))
## merge condition 4
dataset_merged = merge(dataset_merged, dataset[dataset$ms_internal_id==str_join("PTM_UBI_",target,"_IFN_MG132"),], by="peptide", all.x=T, all.y=T)
dataset_merged = dataset_merged[, c(1:7,10,11)] 
colnames(dataset_merged)[8:9] = c(str_join("GFP_H_",target,"_L_IFN_MG132"),str_join(target,"_H_GFP_L_IFN_MG132"))
## chop off first column to prepare for limma and put it as rowname
rownames(dataset_merged) = dataset_merged[,1]
dataset_merged = dataset_merged[,2:9]


###########################
## QUALITY CONTROL

colnames(dataset_merged) = lapply(colnames(dataset_merged),function (x) str_join(substr(x,1,11),"\n",substr(x,13,30)))
setwd(plotDir)
pdf(str_join(target,"_quality_correlations.pdf"), width=8,height=8)
pairs(dataset_merged,upper.panel=panel.cor, main = "", lower.panel=function(x,y)points(x,y,pch=20))  
dev.off()

###########################
## LIMMA

#####################
##  do not modify here

doLimma = function(filename, modeldata, design, biolrep){
  #corfit <- duplicateCorrelation(modeldata, design, ndups = 1, block = biolrep)
  fit <- lmFit(modeldata, design)
  fit <- eBayes(fit)
  result = topTable(fit, adjust = "BH", sort.by="p", number=1000000)
  result = cbind(1:nrow(result),result)
  colnames(result)[1] = "hit_rank"
  result
}

###############################################################
## select some specific sets here based on quality control plots
## and adjust design and biolrep definition
## leave out NA's or not

setwd(processedDataDir)
setwd(str_join(processedDataDir,"/LIMMA_NORMALIZED"))

## write loop for conditions
conditions = c("UNTREATED", "MG132", "IFN", "IFN_MG132")
i = 1
for(c in conditions){
  filename = str_join(target,"_",c,"_limma.txt")
  print(filename)
  modeldata = dataset_merged[,((i*2)-1):(i*2)] 
  modeldata = na.omit(modeldata)
  biolrep = c(1,2)
  design <- c(1,-1)
  result = doLimma(filename, modeldata, design, biolrep)
  # join with peptide-protein lookup table to get readable output and detect if multiple peptides from the same protein were detected
  result_annotated = sqldf(str_join("select * from result R left join peptide_protein_descriptions P on P.peptide = R.ID group by P.peptide, P.key order by P_Value asc"))
  write.table(result_annotated,file=filename, sep="\t", quote=F, eol="\n",row.names=F, col.names=T)
  i = i+1
}

## VIF ###########################################

# filename = "VIF_untreated_limma.txt"
# modeldata = dataset_merged[,1:2] 
# modeldata = na.omit(modeldata)
# biolrep = c(1,2)
# design <- c(1,-1)
# result = doLimma(filename, modeldata, design, biolrep)
# # join with peptide-protein lookup table to get readable output and detect if multiple peptides from the same protein were detected
# result_annotated = sqldf(str_join("select * from result R left join peptide_protein_descriptions P on P.peptide = R.ID group by P.peptide, P.key order by P_Value asc"))
# result_annotated = getUniprotDescriptions(result_annotated, key="proteins", splitkeychar=';')
# write.table(result_annotated,file=filename, sep="\t", quote=F, eol="\n",row.names=F, col.names=T)

