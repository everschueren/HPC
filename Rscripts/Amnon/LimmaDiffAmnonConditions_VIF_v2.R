library(sqldf)
library(plyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(limma)

setwd(scriptDir)
source("Plotting/Aux.R")

PVALUE=0.05

#############################################################################################################################
## 1. COMPARING CONDITIONS FOR VIF

#############################################################################################################################
## 1.1 USING NR> OF IDENTIFIED PEPTIDES

input <- read.delim("~/Projects/HPCKrogan/Amnon/input/vif_rep2/041213-ag-13-24-peptide-counts-reorg.txt")
colnames(input)[1] = "ID"

datamatrix = input[,3:14]
rownames(datamatrix) = input[,1]
colnames(datamatrix) = c("T1_INFECTED","T2_INFECTED","T3_INFECTED","T4_INFECTED","T1_CONTROL","T2_CONTROL","T3_CONTROL","T4_CONTROL","T1_UNINFECTED","T2_UNINFECTED","T3_UNINFECTED","T4_UNINFECTED")

design = model.matrix(~ 0 + factor(rep(1:3,each=4)))
colnames(design) = c("INFECTED", "CONTROL", "UNINFECTED")
fit = lmFit(datamatrix,design)

contrasts = makeContrasts(INFECTED-CONTROL,INFECTED-UNINFECTED, levels=design)
contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
testResults = decideTests(contrasts.fit, method="global")
summary(testResults)

diff_infected_control = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="P", p.value=PVALUE, number=Inf)
diff_infected_control_annotated = merge(diff_infected_control, input[,1:3],by="ID")
diff_infected_control_annotated = sqldf("select ID, logFC, adj_P_Val, B, Fasta_headers from diff_infected_control_annotated where logFC>0 order by logFC desc")

diff_infected_uninfected = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", p.value=PVALUE, number=Inf)
diff_infected_uninfected_annotated = merge(diff_infected_uninfected, input[,1:3],by="ID")
diff_infected_uninfected_annotated = sqldf("select ID, Proteins, logFC, adj_P_Val, B, Fasta_headers from diff_infected_uninfected_annotated where logFC>0 order by logFC desc")

contrasts = makeContrasts(INFECTED-(CONTROL+UNINFECTED), levels=design)
contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
testResults = decideTests(contrasts.fit, method="global")
summary(testResults)
diff = topTable(contrasts.fit, adjust.method="BH", p.value=PVALUE, number=Inf)
diff_annotated = merge(diff, input[,1:3],by="ID")
diff_annotated = sqldf("select ID, logFC, adj_P_Val, B, Fasta_headers from diff_annotated where logFC>0 order by logFC desc")

#############################################################################################################################
## 1.2 USING INTENSITIES

input <- read.delim("~/Projects/HPCKrogan/Amnon/input/vif_rep2/041213-ag-13-24-peptides-intensity-reorg.txt")
colnames(input)[1] = "ID"
datamatrix = input[,4:15]
rownames(datamatrix) = input[,1]
colnames(datamatrix) = c("T1_INFECTED","T2_INFECTED","T3_INFECTED","T4_INFECTED","T1_CONTROL","T2_CONTROL","T3_CONTROL","T4_CONTROL","T1_UNINFECTED","T2_UNINFECTED","T3_UNINFECTED","T4_UNINFECTED")
datamatrix[datamatrix==0] = 1
datamatrix = log2(datamatrix)

## plot correlation
for(i in 1:3){
  pairs(datamatrix[,rep(1:4)+((i-1)*4)], upper.panel=panel.cor)  
}

design = model.matrix(~ 0 + factor(rep(1:3,each=4)))
colnames(design) = c("INFECTED", "CONTROL", "UNINFECTED")
biolrep = rep(1:3,each=4)

fit = lmFit(datamatrix,design)

contrasts = makeContrasts(INFECTED-CONTROL,INFECTED-UNINFECTED, levels=design)
contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
testResults = decideTests(contrasts.fit, method="global")
summary(testResults)

diff_infected_control = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="P", p.value=PVALUE, number=Inf)
diff_infected_control_annotated = merge(diff_infected_control, input[,1:3],by="ID")
diff_infected_control_annotated = sqldf("select ID, Proteins, logFC, adj_P_Val, B, Fasta_headers from diff_infected_control_annotated where logFC>0 order by logFC desc")

diff_infected_uninfected = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", p.value=PVALUE, number=Inf)
diff_infected_uninfected_annotated = merge(diff_infected_uninfected, input[,1:3],by="ID")
diff_infected_uninfected_annotated = sqldf("select ID, logFC, adj_P_Val, B, Fasta_headers from diff_infected_uninfected_annotated where logFC>0 order by logFC desc")

contrasts = makeContrasts(INFECTED-(CONTROL+UNINFECTED), levels=design)
contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
testResults = decideTests(contrasts.fit, method="global")
summary(testResults)
diff = topTable(contrasts.fit, adjust.method="BH", p.value=PVALUE, number=Inf)
diff_annotated = merge(diff, input[,1:3],by="ID")
diff_annotated = sqldf("select ID, Proteins, logFC, adj_P_Val, B, Fasta_headers from diff_annotated where logFC>0 order by logFC desc")




