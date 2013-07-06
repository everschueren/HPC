library(sqldf)
library(plyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(limma)

scriptDir = "/Users/everschueren/Projects/HPCKrogan/Scripts/Rscripts/"
plotDir = "/Users/everschueren/Projects/HPCKrogan/Amnon/output/"

setwd(scriptDir)
source("Plotting/Aux.R")
source("DB/LoadProspector.R")


prospectorData = LoadProspectorData()

###########################################################
## estimate FPR

decoys = sqldf("select ms_machine, count(*) as 'count' from prospectorData where uniprot_ac ='decoy' group by ms_machine")
all = sqldf("select ms_machine, count(*) as 'count' from prospectorData group by ms_machine")

barplot(decoys[,2]/all[,2], ylim=c(0,.1),names.arg=all[,1], main="FDR between Calcium and Elite", ylab= "FDR (%)")

prospectorSelected = sqldf("select ms_entry_id, ms_machine, uniprot_id, unique_peptides_num, cell_timepoint, cell_treatment from prospectorData where uniprot_ac != 'decoy'")

###########################################################
## Large category exploration
prospectorNumIdentified = sqldf("select ms_machine, cell_timepoint, cell_treatment, count(*) as 'count', avg(unique_peptides_num) as mean_unique_peptides_num, std(unique_peptides_num) as se_unique_peptides_num from prospectorSelected group by ms_machine, cell_timepoint, cell_treatment")

# boxplot(count ~ cell_treatment, data=prospectorNumIdentified[prospectorNumIdentified$ms_machine=="Calcium",])
# boxplot(count ~ cell_treatment, data=prospectorNumIdentified[prospectorNumIdentified$ms_machine=="Elite",], add=T)

setwd(plotDir)
pdf("machine_comparison_conditions_countidentified.pdf", width=5, height=5)
qplot(factor(cell_treatment),data=prospectorNumIdentified,geom="bar",fill=ms_machine,weight=count,position="dodge", xlab="",ylab="#Identified proteins", main="Machine comparison over samples") + theme(axis.text.x=element_text(angle=90))
dev.off()

pdf("machine_comparison_timepoints_countidentified.pdf", width=5, height=5)
qplot(factor(cell_timepoint),data=prospectorNumIdentified,geom="bar",fill=ms_machine,weight=count,position="dodge", xlab="",ylab="#Identified proteins", main="Machine comparison over timepoints") 
dev.off()

pdf("machine_comparison_samples_timesidentified_means.pdf", width=5, height=5)
qplot(factor(cell_treatment),data=prospectorNumIdentified,geom="bar",fill=ms_machine,weight=mean_unique_peptides_num,position="dodge", xlab="",ylab="Mean abundance Identified proteins", main="Machine comparison over samples") + theme(axis.text.x=element_text(angle=90))
dev.off()

pdf("machine_comparison_timepoints_timesidentified_means.pdf", width=5, height=5)
qplot(factor(cell_timepoint),data=prospectorNumIdentified,geom="bar",fill=ms_machine,weight=mean_unique_peptides_num,position="dodge", xlab="",ylab="Mean abundance Identified proteins", main="Machine comparison over timepoints")
dev.off()


###########################################################
## correlation between all repeated experiments on two machines

machine_replicas = sqldf("select * from prospectorSelected R1 join prospectorSelected R2 on R1.uniprot_id = R2.uniprot_id where R1.ms_machine = 'Elite' and R2.ms_machine = 'Calcium' and R1.cell_timepoint = R2.cell_timepoint and R1.cell_treatment = R2.cell_treatment")
colnames(machine_replicas)[7:12] = str_join(colnames(machine_replicas)[1:6],"_r")
machine_replicas = sqldf("select * from machine_replicas order by unique_peptides_num desc")

machine_correlation = cor(log2(machine_replicas$unique_peptides_num),log2(machine_replicas$unique_peptides_num_r),method="pearson")
pdf("machine_comparison_overall.pdf", width=5, height=5)
plot(log2(machine_replicas$unique_peptides_num), log2(machine_replicas$unique_peptides_num_r), ylab="log2-abundance (Calcium)", xlab="log2-abundance (Elite)", main=str_join("Correlation of #Identified peptides/protein = ",round(machine_correlation,2)), pch=20)
dev.off()

###########################################################
## correlation between time-points under a single condition
## only shared proteins between conditions

# ## different conditions
# condition.hiv.calcium.q = "select * from prospectorSelected where cell_treatment = 'hiv_infected' and ms_machine='calcium'"
# condition.hivD.calcium.q = "select * from prospectorSelected where cell_treatment = 'HIV_D_vif_infected' and ms_machine='calcium'"
# condition.untreated.calcium.q = "select * from prospectorSelected where cell_treatment = 'untreated' and ms_machine='calcium'"
# condition.hiv.elite.q = "select * from prospectorSelected where cell_treatment = 'hiv_infected' and ms_machine='elite'"
# condition.hivD.elite.q = "select * from prospectorSelected where cell_treatment = 'HIV_D_vif_infected' and ms_machine='elite'"
# condition.untreated.elite.q = "select * from prospectorSelected where cell_treatment = 'untreated' and ms_machine='elite'"
# 
# conditions = c(condition.hiv.calcium.q,condition.hivD.calcium.q,condition.untreated.calcium.q,condition.hiv.elite.q,condition.hivD.elite.q,condition.untreated.elite.q)
# 
# ## select condition data, join timepoints and plot
# for (c in conditions){
#   sample = sqldf(c)
#   sample.timepoints = sqldf("select T1.uniprot_id, T1.unique_peptides_num as 'T1', T2.unique_peptides_num as 'T2', T3.unique_peptides_num as 'T3', T4.unique_peptides_num as 'T4' from sample T1 left join sample T2 on T1.uniprot_id = T2.uniprot_id left join sample T3 on T1.uniprot_id = T3.uniprot_id left join sample T4 on T1.uniprot_id = T4.uniprot_id where T1.cell_timepoint = 6 and T2.cell_timepoint = 16 and T3.cell_timepoint = 26 and T4.cell_timepoint = 42")
#   pairs(sample.timepoints[,2:5],upper.panel=panel.cor, main=substr(c, 40, nchar(c)))  
# }

###########################################################
## correlation between conditions under a single time-point
## only shared proteins between conditions
# 
# condition.T1.calcium.q = "select * from prospectorSelected where cell_timepoint = 6 and ms_machine='calcium'"
# condition.T2.calcium.q = "select * from prospectorSelected where cell_timepoint = 16 and ms_machine='calcium'"
# condition.T1.elite.q = "select * from prospectorSelected where cell_timepoint = 6 and ms_machine='elite'"
# condition.T2.elite.q = "select * from prospectorSelected where cell_timepoint = 16 and ms_machine='elite'"
# 
# ## select condition data, join timepoints and plot
# conditions = c(condition.T1.calcium.q, condition.T2.calcium.q, condition.T1.elite.q, condition.T2.elite.q)
# for (c in conditions){
#   sample = sqldf(c)
#   sample.conditions = sqldf("select S1.uniprot_id, S1.unique_peptides_num as 'S1', S2.unique_peptides_num as 'S2', S3.unique_peptides_num as 'S3' from sample S1 left join sample S2 on S1.uniprot_id = S2.uniprot_id left join sample S3 on S1.uniprot_id = S3.uniprot_id where S1.cell_treatment = 'untreated' and S2.cell_treatment = 'HIV_infected' and S3.cell_treatment ='HIV_D_vif_infected'")
#   pairs(sample.conditions[,2:4],upper.panel=panel.cor, main=substr(c, 40, nchar(c)))
# }

###########################################################
## correlation between conditions under a single time-point and time-points under a single condition
## all proteins, incluidng those not shared between conditions

machines = c("elite","calcium")
timepoints = c(6,16,26,42)
conditions = c("uninfected","HIV_infected","HIV_infected_control")

# LOOP OVER TIMEPOINTS AND COMPARE CONDITIONS

par(pch=20)
for(m in machines){
  for( t in timepoints){
    q = str_join("select * from prospectorSelected where cell_timepoint = ",t," and ms_machine='",m,"'")
    sample = as.data.frame(sqldf(q))
    c1 = sample[sample$cell_treatment=="uninfected",]
    c2 = sample[sample$cell_treatment=="HIV_infected",]
    c3 = sample[sample$cell_treatment=="HIV_infected_control",]
    conditions_treatments = merge(c1, c2, by="uniprot_id", all.x=T, all.y=T)
    conditions_treatments = merge(conditions_treatments, c3, by="uniprot_id", all.x=T, all.y=T)
    
    ## replace NA's by 0 (or 1 for log-scale)
    
    conditions_treatments = conditions_treatments[,c(1,5,4,9,14)]
    conditions_treatments[ is.na(conditions_treatments) ] = 1
    conditions_treatments[,3:5] = log2(conditions_treatments[,3:5])
    
    colnames(conditions_treatments) = c("uniprot_id","timepoint","Uninfected","Infected","Infected-Control")
    overlaps = vennCounts(conditions_treatments[,3:5])
    
    pdf(str_join(m,"_T_",t,"_pairplot.pdf"),width=6,height=6)
    pairs(conditions_treatments[,3:5],upper.panel=panel.cor, main = str_join("Machine=",m," Timepoint=",t),  lower.panel=function(x,y)points(x,y,pch=20))  
    dev.off()
    
    pdf(str_join(m,"_T_",t,"_overlap.pdf"),width=6,height=6)
    vennDiagram(overlaps, main=str_join("Machine=",m," Timepoint=",t))
    dev.off()
  }
}

# LOOP OVER CONDITIONS AND COMPARE TIMEPOINTS
for(m in machines){
  for( c in conditions){
    q = str_join("select * from prospectorSelected where cell_treatment='",c,"' and ms_machine='",m,"'")
    sample = as.data.frame(sqldf(q))
    c1 = sample[sample$cell_timepoint==6,]
    c2 = sample[sample$cell_timepoint==16,]
    c3 = sample[sample$cell_timepoint==26,]
    c4 = sample[sample$cell_timepoint==42,]
    conditions_treatments = merge(c1, c2, by="uniprot_id", all.x=T, all.y=T)
    conditions_treatments = conditions_treatments[,c(1,5,4,9)]
    conditions_treatments = merge(conditions_treatments, c3, by="uniprot_id", all.x=T, all.y=T)
    conditions_treatments = conditions_treatments[,c(1,2,3,4,7)]
    colnames(conditions_treatments) = c("uniprot_id","timepoint","6hrs","16hrs","26hrs")
    conditions_treatments = merge(conditions_treatments, c4, by="uniprot_id", all.x=T, all.y=T)
    conditions_treatments = conditions_treatments[,c(1,2,3,4,5,8)]
    colnames(conditions_treatments) = c("uniprot_id","timepoint","6hrs","16hrs","26hrs","42hrs")
    
    ## replace NA's by 0 (or 1 for log-scale)
    #conditions_treatments[ is.na(conditions_treatments) ] = 0
    conditions_treatments[ is.na(conditions_treatments) ] = 1
    conditions_treatments[,3:6] = log2(conditions_treatments[,3:6])
  
    # overlap
    overlaps = vennCounts(conditions_treatments[,3:5])
    
    pdf(str_join(m,"_C_",c,"_pairplot.pdf"),width=6,height=6)
    pairs(conditions_treatments[,3:6],upper.panel=panel.cor, main = str_join("Machine=",m," Condition=",c), lower.panel=function(x,y)points(x,y,pch=20))  
    dev.off()
    
    pdf(str_join(m,"_C_",c,"_overlap.pdf"),width=6,height=6)
    vennDiagram(overlaps, main = str_join("Machine=",m," Condition=",c))
    dev.off()
  } 
}

