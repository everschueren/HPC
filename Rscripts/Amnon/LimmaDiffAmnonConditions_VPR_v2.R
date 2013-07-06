library(sqldf)
library(plyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(limma)

source("Plotting/Aux.R")
uniprot = read.delim("/Users/everschueren/Projects/HPCKrogan/Data/Uniprot/uniprot_protein_descriptions.txt", stringsAsFactors=F)
jager = read.delim("/Users/everschueren/Projects/HPCKrogan/Data/Jager/nature10719-s3_jurkat.txt", stringsAsFactors=F)  

PVALUE=0.001

#############################################################################################################################
## 1. COMPARING CONDITIONS FOR VPR

#############################################################################################################################
## 1.1 USING NR> OF IDENTIFIED PEPTIDES

input <- read.delim("~/Projects/HPCKrogan/Amnon/input/041313-ag-25-36-peptide-counts-reorg.txt")
colnames(input)[1] = "ID"

datamatrix = input[,3:14]
rownames(datamatrix) = input[,1]
colnames(datamatrix) = c("T1_INFECTED","T2_INFECTED","T3_INFECTED","T4_INFECTED","T1_CONTROL","T2_CONTROL","T3_CONTROL","T4_CONTROL","T1_UNINFECTED","T2_UNINFECTED","T3_UNINFECTED","T4_UNINFECTED")

## plot correlation
for(i in 1:3){
  pairs(datamatrix[,rep(1:4)+((i-1)*4)], upper.panel=panel.cor)  
}

design = model.matrix(~ 0 + factor(rep(1:3,each=4)))
colnames(design) = c("INFECTED", "CONTROL", "UNINFECTED")
fit = lmFit(datamatrix,design)

contrasts = makeContrasts(INFECTED-CONTROL,INFECTED-UNINFECTED, (INFECTED-CONTROL)+(INFECTED-UNINFECTED),levels=design)
contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
testResults = decideTests(contrasts.fit, method="global")
summary(testResults)

diff_infected_control = topTable(contrasts.fit, coef=1, adjust.method="BH", sort.by="P", p.value=PVALUE, number=Inf)
diff_infected_control_annotated = merge(diff_infected_control, input[,1:3],by="ID")
diff_infected_control_annotated = sqldf("select ID, logFC, adj_P_Val, B, Fasta_headers from diff_infected_control_annotated where logFC>0 order by logFC desc")

diff_infected_uninfected = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", p.value=PVALUE, number=Inf)
diff_infected_uninfected_annotated = merge(diff_infected_uninfected, input[,1:3],by="ID")
diff_infected_uninfected_annotated = sqldf("select ID, logFC, adj_P_Val, B, Fasta_headers from diff_infected_uninfected_annotated where logFC>0 order by logFC desc")

diff = topTable(contrasts.fit, coef=3, adjust.method="BH", p.value=PVALUE, number=Inf)
diff_annotated = merge(diff, input[,1:3],by="ID")
diff_annotated = sqldf("select ID, logFC, adj_P_Val, B, Fasta_headers from diff_annotated where logFC>0 order by logFC desc")

#############################################################################################################################
## 1.2 USING INTENSITIES

input <- read.delim("~/Projects/HPCKrogan/Amnon/input/041313-ag-25-36-peptides-intensity-reorg.txt")
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
diff_infected_control_annotated = sqldf("select ID, Proteins, logFC, adj_P_Val, B from diff_infected_control_annotated where logFC>0 order by logFC desc")

diff_infected_uninfected = topTable(contrasts.fit, coef=2, adjust.method="BH", sort.by="P", p.value=PVALUE, number=Inf)
diff_infected_uninfected_annotated = merge(diff_infected_uninfected, input[,1:3],by="ID")
diff_infected_uninfected_annotated = sqldf("select ID, logFC, adj_P_Val, B, from diff_infected_uninfected_annotated where logFC>0 order by logFC desc")

contrasts = makeContrasts(INFECTED-(CONTROL+UNINFECTED), levels=design)
contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
testResults = decideTests(contrasts.fit, method="global")
summary(testResults)
diff = topTable(contrasts.fit, adjust.method="BH", p.value=PVALUE, number=Inf)
diff_annotated = merge(diff, input[,1:3],by="ID")
diff_annotated = sqldf("select ID, Proteins, logFC, adj_P_Val, B from diff_annotated where logFC>0 order by logFC desc")

################################################################################################################################
## 2. COMPARING TIMEPOINTS FOR VPR

#############################################################################################################################
## 2.2 USING INTENSITIES
# 
# input <- read.delim("~/Projects/HPCKrogan/Amnon/input/041313-ag-25-36-peptides-intensity-reorg.txt")
# colnames(input)[1] = "ID"
# datamatrix = input[,4:15]
# rownames(datamatrix) = input[,1]
# colnames(datamatrix) = c("T1_INFECTED","T2_INFECTED","T3_INFECTED","T4_INFECTED","T1_CONTROL","T2_CONTROL","T3_CONTROL","T4_CONTROL","T1_UNINFECTED","T2_UNINFECTED","T3_UNINFECTED","T4_UNINFECTED")
# datamatrix[datamatrix==0] = 1
# datamatrix = log2(datamatrix)
# 
# ## SPECIFICALLY COMPARE TIMEPOINT 2
# 
# design = model.matrix(~ 0 + factor(c(1,2,1,1)))
# colnames(design) = c("T_I","T2_I")
# biolrep = rep(1:12)
# 
# fit = lmFit(datamatrix[,1:4],design)
# contrasts = makeContrasts(T2_I-T_I, levels=design)
# contrasts.fit = eBayes(contrasts.fit(fit, contrasts))
# testResults = decideTests(contrasts.fit, method="global")
# summary(testResults)
# 
# diff_T2 = topTable(contrasts.fit, adjust.method="BH", number=Inf, p.value = 0.001)
# diff_T2_annotated = merge(diff_T2, input[,1:3],by="ID")
# diff_T2_annotated = sqldf("select ID, Proteins, logFC, adj_P_Val, B, Fasta_headers from diff_T2_annotated where logFC>0 order by logFC desc")
# summ = sqldf("select proteins, avg(logFC), avg(-log(adj_P_Val)) as 'p_val' ,count(*) as 'count' from diff_T2_annotated group by proteins")
# summ2 = sqldf("select *, (S.count+0.0)/U.length as 'relative_count' from uniprot U join summ S on S.Proteins = U.uniprot_id  order by relative_count desc")


################################################################################################################################
## 3. time series

significant_set = diff_infected_control_annotated

tmp  = cbind(rownames(datamatrix), datamatrix[,1:4]) 
colnames(tmp)[1] = "ID"
test = merge(tmp, significant_set, by ="ID")
test = sqldf("select Proteins,avg(T1_INFECTED) as 'T1',avg(T2_INFECTED) as 'T2',avg(T3_INFECTED) as 'T3',avg(T4_INFECTED) as 'T4', count(*) as 'count' from test group by Proteins having count > 1")

actin_control = sqldf("select 'P68133' as 'Proteins', avg(T.T1_INFECTED) as 'T1',avg(T.T2_INFECTED) as 'T2',avg(T.T3_INFECTED) as 'T3',avg(T.T4_INFECTED) as 'T4', count(*) as 'count' from tmp T join input I on I.ID = T.ID where I.Proteins like '%P68133%' group by I.proteins having count > 10")

general_control = sqldf("select I.Proteins, avg(T.T1_INFECTED) as 'T1',avg(T.T2_INFECTED) as 'T2',avg(T.T3_INFECTED) as 'T3',avg(T.T4_INFECTED) as 'T4', count(*) as 'count' from tmp T join input I on I.ID = T.ID group by I.Proteins having count > 10")
general_control = sqldf("select 'P68133' as 'Proteins', avg(T1) as 'T1',avg(T2) as 'T2',avg(T3) as 'T3',avg(T4) as 'T4', count(*) as 'count' from general_control")

test =rbind(test,general_control)
test = sqldf("select * from test T join uniprot U on U.uniprot_id = T.Proteins")
test$uniprot_ac[nrow(test)] = "BACKGROUND" 

minval = c(min(test[,2:5]))
maxval = c(max(test[,2:5]))

par(mar=c(8,6,6,4))
plot(NULL,xlim=c(1,4),ylim=c(minval,maxval),xlab="Timepoints",ylab="LOG-intensities",axes=F, main="Stat. Significant interactors for VPR over timepoints\n(INFECTED-NEG.CONTROL,P<.001,#peptides>1)",cex.lab=1.3)
cols = rainbow(nrow(test))

for(r in 1:nrow(test)){
  d = as.numeric(test[r,2:5])
  lines(d,col=cols[r], lwd=3)
}
axis(1,at=1:4,labels=c("T1","T2","T3","T4"))
axis(2,at=round(minval):round(maxval),labels=round(minval):round(maxval),las=2)
par(xpd=F)

legend("bottom", inset=c(-2,-1),test$uniprot_ac,fill=cols)