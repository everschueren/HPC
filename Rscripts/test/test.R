library(sqldf)

# input <- read.delim("~/Projects/HPCKrogan/Amnon/processed/041213-ag-13-24-differential.txt", dec=",")
# 
# test1 = sqldf("select uniprot_ac, avg(infected_neg_control_logFC) as 'avg_logFC', count(*) from input where infected_neg_control_adjPVal<0.05 group by uniprot_ac order by avg_logFC desc")
# 
# test2 = sqldf("select uniprot_ac, avg(infected_uninfected_logFC) as 'avg_logFC', count(*) from input where infected_uninfected_adjPVal<0.05 group by uniprot_ac order by avg_logFC desc")
# 
# test3 = sqldf("select * from test1 T1 join test2 T2 on T1.uniprot_ac = T2.uniprot_ac where T1.avg_logFC > 0 and T2.avg_logFC>0")
# 
# 
# scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
# setwd(scriptDir)
# source("Plotting/Aux.R")
# input <- read.delim("~/Projects/HPCKrogan/Data/Amnon/processed/041213-ag-13-24-matrix-max.txt")
# colnames(input)[1] = "ID"

datamatrix = log2(input[,3:14]+1)

rownames(datamatrix) = input[,1]
colnames(datamatrix) = c("T1_INFECTED","T2_INFECTED","T3_INFECTED","T4_INFECTED","T1_CONTROL","T2_CONTROL","T3_CONTROL","T4_CONTROL","T1_UNINFECTED","T2_UNINFECTED","T3_UNINFECTED","T4_UNINFECTED")

## plot correlation
for(i in 1:3){
  pairs(datamatrix[,rep(1:4)+((i-1)*4)], upper.panel=panel.cor, main="Jurkat-Infection-PPI: Vif Timepoints")  
}


# jager_jurkat <- read.delim("~/Projects/HPCKrogan/Data/Jager/nature10719-s3_jurkat.txt")
# 
# jager_jurkat_vif = sqldf("select * from jager_jurkat where Bait='VIF'")
# jager_jurkat_vif_significant = sqldf("select * from jager_jurkat where mist_score > 0.75 and Bait='VIF' order by mist_score desc")
# overlap = sqldf("select * from jager_jurkat_vif J join test_agg T on T.Uniprot_ac__1 = J.Prey where Bait='VIF' order by infected_neg_control_logFC_avg desc")


# test_agg = aggregate(. ~ uniprot_ac, data = test[,c(1,6,8,10)], mean)
# test_agg = sqldf("select uniprot_ac, avg(infected_neg_control_logFC) as 'avg1', avg(-log(infected_neg_control_adjPVal)) as 'avg2', avg(infected_uninfected_logFC) as 'avg3', avg(-log(infected_uninfected_adjPVal)) as 'avg4', count(*) as 'count' from test group by uniprot_ac order by count desc")