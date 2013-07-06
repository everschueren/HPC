## testing fpr the effect of abundance

library(sqldf)
LIMMA_Vpu_IFN_MG132 <- read.delim("~/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/processed/LIMMA_MAXQ/LIMMA_Vpu_IFN_MG132.txt")

res = sqldf("select * from LIMMA_Vpu_IFN_MG132 group by rank")
res = sqldf("select *, count(*) as 'count' from LIMMA_Vpu_IFN_MG132 group by key")
res = cbind(res, res$count/res$length)
colnames(res)[ncol(res)]="numpepperprotein"
res=cbind(res, res$P_Value*(1-res$numpepperprotein)) 
colnames(res)[ncol(res)]="pvalue_times_abundance"
res = sqldf("select * from res order by pvalue_times_abundance asc")