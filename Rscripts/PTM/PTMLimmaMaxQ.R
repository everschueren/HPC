library(sqldf)
library(limma)
library(reshape2)
library(stringr)

scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
plotDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/output/"
processedDataDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/LIMMA_MAXQ/"
inputDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/MAXQ_MATRICES/"

processedPeptideDataDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/"
peptide_protein_descriptions <- read.delim(str_join(processedPeptideDataDir,"peptide_protein_descriptions.txt"))

runFiles = read.delim("/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/MAXQ_MATRICES/runfiles.txt")

#############
## FUNCTIONS

doLimma = function(dataMatrix, design, biolrep){
  corfit <- duplicateCorrelation(dataMatrix, ndups = 1, block = biolrep)
  print(corfit$consensus)
  fit <- lmFit(dataMatrix, design, block = biolrep, cor = corfit$consensus)
  #fit <- lmFit(test, design)
  cont.matrix <- makeContrasts("(hiv_h+hiv_l)-(gfp_h+gfp_l)",levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  result <- decideTests(fit2,method="global")
  print(summary(result))
  result = topTable(fit2, adjust = "BH",number=100000)  
  result
}

#############
## READ INPUT 
#runFiles = runFiles[1,]

for(i in 1:nrow(runFiles)){
  
  ## get input and conditions for exp and read data
  runFile = runFiles$evidence_file[i]
  target = runFiles$target[i]
  condition = runFiles$condition[i]
  print(str_join(target,"_",condition))
  designColumns = str_split(runFiles$design_columns[i],pattern=',')
  dataMatrix = read.delim(str_join(inputDir, runFile))
  
  ## kick out all non-ubiquitylated peptides by looking for (gl) (GlyGly-antibody)
  dataMatrix = sqldf("select * from dataMatrix where mod_seq like '%(gl)%'")
  
  ## put ids as rownames and take log-intensities
  tmp = dataMatrix[, 2:9]
  rownames(tmp) = dataMatrix[, 1]
  dataMatrix = tmp
  dataMatrix = log2(dataMatrix)
  
  ## Normalize
  # dataMatrix = normalizeBetweenArrays(as.matrix(dataMatrix), method="cyclicloess")
  
  ## make design matrix and call Limma wrapper
  design <- model.matrix(~ 0+factor(c(1,1,2,2,3,3,4,4)))
  biolrep = c(1,1,2,2,3,3,4,4)
  colnames(design) = as.vector(designColumns)[[1]]
  result = doLimma(dataMatrix,design,biolrep)
  
  ## kick out all non-ubiquitylated peptides by looking for (gl) (GlyGly-antibody)
  #result = sqldf("select * from result where ID like '%(gl)%'")
  
  ## format for output and write to file
  result = sqldf("select * from result where logFC > 0")
  result = cbind(1:nrow(result),result)
  colnames(result)[1] = "rank"
  result_annotated = sqldf(str_join("select * from result R left join peptide_protein_descriptions P on P.peptide = R.ID group by P.peptide, P.key order by R.P_Value asc"))
  write.table(result_annotated,file=str_join(processedDataDir,"LIMMA_",target,"_",condition,".txt"), eol="\n", sep="\t", quote=F, row.names=F)
}

