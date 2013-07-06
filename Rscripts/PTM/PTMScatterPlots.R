library(ggplot2)
library(sqldf)
library(stringr)
library(limma)

scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
plotDir = "/Users/everschueren/Projects/HPCKrogan/PTM/HIVAccesoryUbiScan/output/"

setwd(scriptDir)
source("DB/LoadMaxQuant.R")

MaxQData = LoadMaxQuantData()

#####
# OLD

targets = c(20261, 20262, 20263)  #VIF, VPR, VPU
conditions = c("untreated", "MG132", "IFN", "IFN;MG132")

setwd(plotDir)
for(t in targets){
  pdf(str_join(t,"_SILAC_ratios.pdf"),width=8, height=9)
  par(mfrow=c(2,2))
  for (c in conditions){
    targetRatios = sqldf(str_join("select proteins, peptide, sample, silac_heavy, silac_light, avg(avg_log_ratio) as 'avg_log_ratio' from MaxQData where peptide like '%(gl)%' and cell_overexpressed_proteins =", t," and cell_treatment='",c,"' group by peptide, sample"))
    temp = merge(targetRatios[targetRatios$silac_heavy == "GFP",], targetRatios[targetRatios$silac_light == "GFP",], by="peptide")
    plot(temp[,6], temp[,11], pch=20, main=c, xlab=temp$sample.x[1], ylab=temp$sample.y[1])
    abline(h=0)
    abline(v=0)
  }
  dev.off()
}