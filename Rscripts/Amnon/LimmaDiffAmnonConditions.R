library(sqldf)
library(plyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(limma)

scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
plotDir = "/Users/everschueren/Projects/HPCKrogan/Amnon/output/"
processedDataDir = "/Users/everschueren/Projects/HPCKrogan/Amnon/processed/"

setwd(scriptDir)
source("Plotting/Aux.R")
source("DB/LoadProspector.R")

prospectorData = LoadProspectorData()
prospectorSelected = sqldf("select ms_entry_id, ms_internal_id, ms_machine, uniprot_id, unique_peptides_num, cell_timepoint, cell_treatment from prospectorData where uniprot_ac != 'decoy'")

## RESHAPE TO MAKE QUICK SELECTIONS

prospectorSelectedValues = dcast(prospectorSelected, uniprot_id ~ ms_internal_id , value.var=c('unique_peptides_num'), margins=T, fill=0,fun=mean)

######### 
### LIMMA

## format the data-matrix for Limma methods
prospectorSelectedValues = dcast(prospectorSelected, uniprot_id ~ ms_internal_id , value.var=c('unique_peptides_num'), margins=T, fill=0,fun=mean)
rownames(prospectorSelectedValues) = prospectorSelectedValues[,1]
prospectorSelectedValues = prospectorSelectedValues[1:(nrow(prospectorSelectedValues)-1),c(2,6,7,8,9,10,11,12,13,3,4,5,14,18,19,20,21,22,23,24,25,15,16,17)]
## log-transform (+ 1 to avoid -Inf after logarithm)
prospectorSelectedValues = log2(prospectorSelectedValues+1)

## strip _HUMAN suffix 
rownames(prospectorSelectedValues) = as.vector(sapply(rownames(prospectorSelectedValues), function(x) strsplit(x, '_')[[1]])[1,])

################################################
## set some variables describing the experiments 

conditions = c("infected", "control", "uninfected")
diff_infected_control_name = "diff_infected_control"
diff_infected_stable_name = "diff_infected_stable"
diff_infected_vs_stable_and_control_name = "diff_infected_vs_stable_and_control"

timepoints = c(6,16,26,42)
PVALUE = 0.01
#evidence = "calcium" ## c("calcium","elite","both")
#evidence = "elite" ## c("calcium","elite","both")
evidence = "both" ## c("calcium","elite","both")

diff_infected_control_name = str_join(diff_infected_control_name, "_", evidence) 
diff_infected_stable_name = str_join(diff_infected_stable_name, "_", evidence)
diff_infected_vs_stable_and_control_name = str_join(diff_infected_vs_stable_and_control_name, "_", evidence)

################################################
# select subsets from data based on #replicates
if(evidence=="both"){           ## use all evidence, machine-repeats and timepoints are all considered as replicas
  dataMatrix = prospectorSelectedValues
  replicates = 2
}else if(evidence=="calcium"){  ## use only i machine's evidence, timepoints are all considered as replicas
  dataMatrix = prospectorSelectedValues[,1:(length(conditions)*length(timepoints))]
  replicates = 1
}else if(evidence=="elite"){
  dataMatrix = prospectorSelectedValues[,((length(conditions)*length(timepoints))+1):((length(conditions)*length(timepoints))*2)]
  replicates = 1
}

# limma part 
design = model.matrix(~ 0+factor(rep(1:length(conditions),length(timepoints)*replicates)))
colnames(design) = conditions

# biolrep = rep(1:3,8)
# corfit <- duplicateCorrelation(dataMatrix,design,block=biolrep)

fit = lmFit(dataMatrix,design)

contrasts = makeContrasts(str_join(conditions[1],"-",conditions[2]),str_join(conditions[1],"-",conditions[3]), str_join(conditions[1],"-(",conditions[2],"+",conditions[3],")"), levels=design)
contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
testResults = decideTests(contrasts.fit, method="global")
summary(testResults)

####################################################
## sort results by p-values and write output
diff_infected_control = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="P", number=5000)
diff_infected_control = diff_infected_control[diff_infected_control$logFC > 1 & diff_infected_control$adj.P.Val < 0.05, ]

diff_infected_stable = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", number=5000)
diff_infected_stable = diff_infected_stable[diff_infected_stable$logFC > 1 & diff_infected_stable$adj.P.Val < 0.05, ]

diff_infected_vs_stable_and_control = topTable(contrasts.fit, coef=3, adjust.method="BH", sort.by="P", number=5000)
diff_infected_vs_stable_and_control = diff_infected_vs_stable_and_control[diff_infected_vs_stable_and_control$logFC > 1 & diff_infected_vs_stable_and_control$adj.P.Val < 0.05, ]

setwd(processedDataDir)
write.table(diff_infected_control, file=str_join(diff_infected_control_name,".txt"),sep="\t",eol="\n",row.names=F,col.names=T, quote=F)
write.table(diff_infected_stable, file=str_join(diff_infected_stable_name,".txt"),sep="\t",eol="\n",row.names=F,col.names=T, quote=F)
write.table(diff_infected_vs_stable_and_control, file=str_join(diff_infected_vs_stable_and_control_name,".txt"),sep="\t",eol="\n",row.names=F,col.names=T, quote=F)


####################################################
## sort results by log-odds values and write output
# diff_infected_control = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="B", number=5000)
# diff_infected_stable = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="B", number=5000)

# 
# ####################################################
# ## sort results by log-fold-changes and write output
# setwd(plotDir)
# diff_infected_control = topTable(contrasts.fit, coef=1, sort.by="logFC", adjust.method="BH", p.value = PVALUE, number=50)
# pdf(str_join(diff_infected_control_name,".pdf"), width=5, height=5)
# par(mar=c(8,5,4,2))
# barplot(diff_infected_control$logFC[diff_infected_control$logFC>0], names.arg = diff_infected_control$ID[diff_infected_control$logFC>0], beside=T,  las=2, main=str_join("HIV-Vif interactors\ninfected / control cells (p < ",PVALUE,")"), ylab="Log-Fold-Change", cex.lab=1.5)
# dev.off()
# 
# diff_infected_stable = topTable(contrasts.fit, coef=2, sort.by="logFC", adjust.method="BH", p.value = PVALUE, number=50)
# pdf(str_join(diff_infected_stable_name,".pdf"), width=5, height=5)
# par(mar=c(8,5,4,2))
# barplot(diff_infected_stable$logFC[diff_infected_stable$logFC>0], names.arg = diff_infected_stable$ID[diff_infected_stable$logFC>0], beside=T,  las=2, main=str_join("HIV-Vif interactors\ninfected / uninfected cells (p < ",PVALUE,")"), ylab="Log-Fold-Change", cex.lab=1.5)
# dev.off()
# 
# 
