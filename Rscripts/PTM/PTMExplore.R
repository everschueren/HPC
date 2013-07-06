library(ggplot2)
library(sqldf)
source("PTMLoad.R")

curdir = getwd()
## use sqldf to query data
overview = sqldf("select R.ms_experiment_id, E.cell_treatment, R.silac_heavy, R.silac_light, E.cell_overexpressed_proteins, count(*) as 'count' from MS_Result_simple R join MS_Experiment E on R.ms_experiment_id=E.id where peptide like '%(gl)%' group by R.ms_experiment_id, R.sample", stringsAsFactors=F)

## OK

ptm_overview_theme = theme_update(axis.text = element_text(size=16), axis.title = element_text(size=16), panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(), legend.text=element_text(size=16), legend.title=element_text(size=16))

qplot(ms_experiment_id, count,data=overview, color=cell_treatment, geom=c("boxplot")) + facet_grid(. ~cell_overexpressed_proteins) + xlab("Overexpressed HIV Proteins") + ylab("# Identified peptides")

setwd("../../PTM/HIVAccesoryUbiScan/output/")
ggsave("Indentified_peptides_across_experiments.pdf")
setwd(curdir)
