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
setwd(processedDataDir)
write.table(prospectorSelectedValues, file="alldata_combined.txt", eol="\n", quote=F, sep="\t", row.names=T, col.names=T)

## log-transform (+ 1 to avoid -Inf after logarithm)
prospectorSelectedValues = log2(prospectorSelectedValues+1)
write.table(total, file="diff_infected_vs_control_and_stable_per_timepoint.txt", eol="\n", quote=F, sep="\t", row.names=F, col.names=F)
## strip _HUMAN suffix 
#rownames(prospectorSelectedValues) = as.vector(sapply(rownames(prospectorSelectedValues), function(x) strsplit(x, '_')[[1]])[1,])



################################################
## set some variables describing the experiments 

timepoints = c(6,16,26,42)
PVALUE = 0.01

# ################################
# ## DATASET : INFECTED VS CONTROL
# 
# timepointlabels =c("infhr06","infhr16","infhr26","infhr42", "infctrlhr06","infctrlhr16","infctrlhr26","infctrlhr42")
# ## calcium 
# inf_ctrl_timepoints = c(1,4,7,10,2,5,8,11)
# ## adds elite repeats
# dataMatrix = prospectorSelectedValues[,c(inf_ctrl_timepoints,inf_ctrl_timepoints+12)]
# design = model.matrix(~0 + factor(rep(1:8,2)))
# colnames(design) = timepointlabels
# # biolrep = rep(1:8,2)
# # corfit <- duplicateCorrelation(dataMatrix,design,block=biolrep)
# fit = lmFit(dataMatrix,design)
# 
# ##############################################################################################################################
# ##############################################################################################################################
# ## contrasts between infected timepoints
# ## => rubbish
# 
# ## diff specific at 1 st timepoint
# contrasts = makeContrasts(str_join(timepointlabels[2],"-",timepointlabels[1]), str_join(timepointlabels[3],"-",timepointlabels[2]), str_join(timepointlabels[3],"-",timepointlabels[1]), str_join(timepointlabels[4],"-",timepointlabels[3]), str_join(timepointlabels[4],"-",timepointlabels[2]), str_join(timepointlabels[4],"-",timepointlabels[1]), levels=design)
# contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
# testResults = decideTests(contrasts.fit, method="global")
# 
# testResultSummary = summary(testResults)
# par(mar=c(8,4,4,2))
# barplot(testResultSummary[c(1,3),], las=2, beside=T)
# 
# diff_infected_16_06 = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="P", number=5000)
# diff_infected_26_16 = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", number=5000)
# diff_infected_26_06 = topTable(contrasts.fit, coef=3, adjust.method="BH", sort.by="P", number=5000)
# diff_infected_42_26 = topTable(contrasts.fit, coef=4, adjust.method="BH", sort.by="P", number=5000)
# diff_infected_42_16 = topTable(contrasts.fit, coef=5, adjust.method="BH", sort.by="P", number=5000)
# diff_infected_42_06 = topTable(contrasts.fit, coef=6, adjust.method="BH", sort.by="P", number=5000)
# 
# ##############################################################################################################################
# ##############################################################################################################################
# ## contrasts between control and infected on separate timepoints
# ## => useful
# 
# contrasts = makeContrasts(str_join(timepointlabels[1],"-",timepointlabels[5]), str_join(timepointlabels[2],"-",timepointlabels[6]), str_join(timepointlabels[3],"-",timepointlabels[7]), str_join(timepointlabels[4],"-",timepointlabels[8]), levels=design)
# contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
# testResults = decideTests(contrasts.fit, method="global")
# 
# testResultSummary = summary(testResults)
# par(mar=c(8,4,4,2))
# barplot(testResultSummary[c(1,3),], las=2, beside=T)
# 
# diff_infected_control_06 = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="P", number=1000)
# diff_infected_control_16 = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", number=1000)
# diff_infected_control_26 = topTable(contrasts.fit, coef=3, adjust.method="BH", sort.by="P", number=1000)
# diff_infected_control_42 = topTable(contrasts.fit, coef=4, adjust.method="BH", sort.by="P", number=1000)
# 
# #############################################
# ## DATASET : INFECTED VS NON-INFECTED (STABLE)
# 
# timepointlabels =c("inf_06","inf_16","inf_26","inf_42", "sta_06","sta_16","sta_26","sta_42")
# ## calcium 
# inf_sta_timepoints = c(1,4,7,10,3,6,9,12)
# dataMatrix = prospectorSelectedValues[,inf_sta_timepoints]
# design = model.matrix(~0 + factor(rep(1:8,1)))
# 
# ## adds elite repeats
# dataMatrix = prospectorSelectedValues[,c(inf_sta_timepoints,inf_sta_timepoints+12)]
# design = model.matrix(~0 + factor(rep(1:8,2)))
# 
# colnames(design) = timepointlabels
# #biolrep = rep(1:8,2)
# #corfit <- duplicateCorrelation(dataMatrix,design,block=biolrep)
# fit = lmFit(dataMatrix,design)
# 
# ##########################################################################################
# ##########################################################################################
# ## contrasts between infected and stable on separate timepoints
# ## => useful
# 
# contrasts = makeContrasts(str_join(timepointlabels[1],"-",timepointlabels[5]), str_join(timepointlabels[2],"-",timepointlabels[6]), str_join(timepointlabels[3],"-",timepointlabels[7]), str_join(timepointlabels[4],"-",timepointlabels[8]), levels=design)
# 
# # contrasts = makeContrasts(str_join(timepointlabels[5],"-",timepointlabels[1]), str_join(timepointlabels[6],"-",timepointlabels[2]), str_join(timepointlabels[7],"-",timepointlabels[3]), str_join(timepointlabels[8],"-",timepointlabels[4]), levels=design)
# contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
# testResults = decideTests(contrasts.fit, method="global")
# 
# testResultSummary = summary(testResults)
# par(mar=c(8,4,4,2))
# barplot(testResultSummary[c(1,3),], las=2, beside=T)
# 
# PVALUE = 0.05
# diff_infected_stable_06 = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="P", p.value =PVALUE, number=1000)
# diff_infected_stable_06 = diff_infected_stable_06[diff_infected_stable_06$logFC>0, ]
# 
# diff_infected_stable_16 = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", p.value =PVALUE)
# diff_infected_stable_16 = diff_infected_stable_16[diff_infected_stable_16$logFC>0, ]
# 
# diff_infected_stable_26 = topTable(contrasts.fit, coef=3, adjust.method="BH", sort.by="P", p.value =PVALUE, number=100)
# diff_infected_stable_42 = topTable(contrasts.fit, coef=4, adjust.method="BH", sort.by="P", p.value =PVALUE, number=100)

#############################################$
## DATASET : INFECTED VS (NON-INFECTED+STABLE)
## best results

timepointlabels =c("inf_06","inf_16","inf_26","inf_42", "ctrl_06","ctrl_16","ctrl_26","ctrl_42", "sta_06","sta_16","sta_26","sta_42")
## calcium 

inf_sta_timepoints = c(1,4,7,10,2,5,8,11,3,6,9,12)

## add elite repeats
dataMatrix = prospectorSelectedValues[,c(inf_sta_timepoints,inf_sta_timepoints+12)]
design = model.matrix(~0 + factor(rep(1:12,2)))
colnames(design) = timepointlabels

#biolrep = rep(1:12,2)
#corfit <- duplicateCorrelation(dataMatrix,design,block=biolrep)
fit = lmFit(dataMatrix,design)

#######################################################################
#######################################################################
## contrasts between infected and stable+control on separate timepoints
## => very useful

contrasts = makeContrasts(str_join(timepointlabels[1],"-(",timepointlabels[5],"+",timepointlabels[9],")"), str_join(timepointlabels[2],"-(",timepointlabels[6],"+",timepointlabels[10],")"), str_join(timepointlabels[3],"-(",timepointlabels[7],"+",timepointlabels[11],")"),str_join(timepointlabels[4],"-(",timepointlabels[8],"+",timepointlabels[12],")"), levels=design)
contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
testResults = decideTests(contrasts.fit, method="global")

testResultSummary = summary(testResults)

PVALUE = 1
LFC = 1
diff_06 = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="P", number=1000)
diff_06 = diff_06[diff_06$logFC>=LFC & diff_06$adj.P.Val<=PVALUE, ]

diff_16 = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", number=1000)
diff_16 = diff_16[diff_16$logFC>=LFC & diff_16$adj.P.Val<=PVALUE, ]

diff_26 = topTable(contrasts.fit, coef=3, adjust.method="BH", sort.by="P", number=1000)
diff_26 = diff_26[diff_26$logFC>=LFC & diff_26$adj.P.Val<=PVALUE, ]

diff_42 = topTable(contrasts.fit, coef=4, adjust.method="BH", sort.by="P", number=1000)
diff_42 = diff_42[diff_42$logFC>=LFC & diff_42$adj.P.Val<=PVALUE, ]

total = merge(diff_06[,c(1,2,6)], diff_16[,c(1,2,6)], by="ID", all.x=T, all.y=T)
total = merge(total, diff_26[,c(1,2,6)], by="ID", all.x=T, all.y=T)
colnames(total) = c("ID","06hr.LFC", "06hr.P","16hr.LFC", "16hr.P","26hr.LFC", "26hr.P")
total = merge(total, diff_42[,c(1,2,6)], by="ID", all.x=T, all.y=T)
colnames(total) = c("ID","06hr.LFC", "06hr.P","16hr.LFC", "16hr.P","26hr.LFC", "26hr.P","42hr.LFC", "42hr.P")
setwd(processedDataDir)
write.table(total, file="diff_infected_vs_control_and_stable_per_timepoint.txt", eol="\n", quote=F, sep="\t", row.names=F, col.names=F)