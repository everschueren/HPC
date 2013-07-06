library(sqldf)
library(limma)
library(reshape2)
library(stringr)
library(mice)

scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
plotDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/output/"
processedDataDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/LIMMA_MAXQ_ALLGFP/"
inputDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/MAXQ_MATRICES/"

processedPeptideDataDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/"
peptide_protein_descriptions <- read.delim(str_join(processedPeptideDataDir,"peptide_protein_descriptions.txt"))

runFiles = read.delim("/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/MAXQ_MATRICES/runfiles.txt", stringsAsFactors=F)

###########################
### COMBINE ALL EXPERIMENTS

technical_replicates = 2
biological_replicates=2
bigMatrix = NULL

for(i in 1:nrow(runFiles)){
  ## get input and conditions for exp and read data
  runFile = runFiles$evidence_file[i]
  target = runFiles$target[i]
  condition = runFiles$condition[i]
  print(str_join(target,"_",condition))
  designColumns = str_split(runFiles$design_columns[i],pattern=',')
  dataMatrix = read.delim(str_join(inputDir, runFile))
  if(is.null(bigMatrix)){
    bigMatrix = dataMatrix
  }else{
    bigMatrix = merge(bigMatrix, dataMatrix, by="mod_seq", all.x=T, all.y=T)
  }
}

matrixColnames = c()
for(r in 1:nrow(runFiles)){
  row = runFiles[r,]
  baseName = str_join(row$target,"_",row$condition)
  designColumns = str_split(row$design_columns, pattern=',')[[1]]
  for(d in 1:length(designColumns)){
    biolName = str_join(baseName,"_", designColumns[d])
    for(t in 1:technical_replicates){
      techName = str_join(biolName, "_", t)
      matrixColnames = c(matrixColnames, techName)
    }
  }
}
colnames(bigMatrix) = c("ID",matrixColnames)

########################################################################
## CONVERT TO LOG VALUES

tmp = bigMatrix[, 2:ncol(bigMatrix)]
rownames(tmp) = bigMatrix[, 1]
bigMatrix = tmp
#bigMatrix[bigMatrix==0 | is.na(bigMatrix)]=1
bigMatrix[bigMatrix==0]=1
bigMatrix = log2(bigMatrix)

########################################################################

#bigMatrix = sqldf("select * from bigMatrix where mod_seq like '%(gl)%'")

## exract all GFP columns
tmp =  as.data.frame(cbind(1:(length(matrixColnames)), matrixColnames))
colnames(tmp) = c("id","name")
gfp_cols = as.numeric(as.character(sqldf("select * from tmp where name like '%gfp%'", stringsAsFactors = FALSE)$id))
gfp_matrix = bigMatrix[, gfp_cols]

write.table(gfp_matrix,file=str_join(processedPeptideDataDir,"GFPLogMatrix.txt"),eol="\n",sep="\t",quote=F, row.names=F, col.names=T)
write.table(gfp_matrix_amelia$imputations[[1]],file=str_join(processedPeptideDataDir,"GFPLogMatrix_amelia.txt"),eol="\n",sep="\t",quote=F, row.names=F, col.names=T)
corMatrix = as.matrix(cor(gfp_matrix, method="spearman"))
heatmap(corMatrix)

gfp_matrix <- read.delim("~/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/GFPLogMatrix_amelia.txt")
gfp_matrix <- read.delim("~/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/GFPLogMatrix.txt")

#############
## FUNCTIONS

doLimma = function(dataMatrix, design, biolrep, contrasts){
  corfit <- duplicateCorrelation(dataMatrix, ndups = 1, block = biolrep)
  print(corfit$consensus)
  fit <- lmFit(dataMatrix, design, block = biolrep, cor = corfit$consensus)
  #fit <- lmFit(test, design)
  
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  result <- decideTests(fit2,method="global")
  print(summary(result))
  result = topTable(fit2, adjust = "BH",number=100000)  
  result
}

## function w. duplicatecorrelation
# bigMatrix = normalizeBetweenArrays(as.matrix(bigMatrix), method="cyclicloess")
# selectMatrix = bigMatrix[,as.vector(c(hivCols,gfpCols))]
# result = doLimma(selectMatrix,design,biolrep)
# result_annotated = sqldf(str_join("select * from result R left join peptide_protein_descriptions P on P.peptide = R.ID group by P.peptide, P.key order by R.B desc"))


# #############
# ## READ INPUT 
#runFiles = runFiles[1,]

for(i in 1:nrow(runFiles)){
  
  ## get input and conditions for exp and read data 
  runFile = runFiles$evidence_file[i]
  target = runFiles$target[i]
  condition = runFiles$condition[i]
  print(str_join(target,"_",condition))
  idx = ((i-1)*8)+1
  dataMatrix = bigMatrix[,c(idx:(idx+7))]
  
  ## filter out hiv cols from bigmatrix
  colIdxs = c()
  for(j in 1:ncol(dataMatrix)){
    colname = colnames(dataMatrix)[j]
    if(grepl("hiv",colname)){
      colIdxs = c(colIdxs, j)
    }
  }
  
  ## filter out gfp cols from gfpmatrix
  colGfpIdxs = c()
  for(j in 1:ncol(gfp_matrix)){
    colname = colnames(gfp_matrix)[j]
    if(grepl(str_join("_",condition,"_"),colname)){
      colGfpIdxs = c(colGfpIdxs, j)
    }
  }
  
  dataWithGfpMatrix = cbind(dataMatrix[,colIdxs], gfp_matrix[,colGfpIdxs])
 # dataWithGfpMatrix = normalizeBetweenArrays(as.matrix(dataWithGfpMatrix), method="quantile")
  design = model.matrix(~ 0+factor(c(rep(1,technical_replicates*biological_replicates),rep(2,ncol(gfp_matrix)/4))))
  colnames(design) = c("hiv","gfp")
  biolrep = rep(1:(ncol(dataWithGfpMatrix)/technical_replicates),each=technical_replicates)
  contrasts <- makeContrasts("hiv-gfp",levels=design)
  result = doLimma(dataWithGfpMatrix,design,biolrep, contrasts)
  result = sqldf("select * from result where logFC > 0")
  result = cbind(1:nrow(result),result)
  colnames(result)[1] = "rank"
  result_annotated = sqldf(str_join("select * from result R left join peptide_protein_descriptions P on P.peptide = R.ID where R.ID like '%gl%' group by P.peptide, P.key order by R.P_Value asc"))
  
  fileName = str_join(processedDataDir,"LIMMA_",target,"_",condition,".txt")
  print(fileName)
  write.table(result_annotated,file=fileName, eol="\n", sep="\t", quote=F, row.names=F)
}

