source("../Rscripts/Plotting/MyPlots.R")

vif_unfiltered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_ratio_intensity/042313-tlj-59-62-matrix-ratios-unfiltered.txt")
vif_ubifiltered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_ratio_intensity/042313-tlj-59-62-matrix-ratios-ubifiltered.txt")
vif_ubi_unique_filtered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_ratio_intensity/042313-tlj-59-62-matrix-ratios-ubi_unique_filtered.txt")

vpr_unfiltered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpr-Ub/processed_ratio_intensity/042313-tlj-57-60-matrix-ratios-unfiltered.txt")
vpr_ubifiltered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpr-Ub/processed_ratio_intensity/042313-tlj-57-60-matrix-ratios-ubifiltered.txt")
vpr_ubi_unique_filtered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpr-Ub/processed_ratio_intensity/042313-tlj-57-60-matrix-ratios-ubi_unique_filtered.txt")

vpu_unfiltered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-Ub/processed_ratio_intensity/051313-tlj-52-55-matrix-ratios-unfiltered.txt")
vpu_ubifiltered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-Ub/processed_ratio_intensity/051313-tlj-52-55-matrix-ratios-ubifiltered.txt")
vpu_ubi_unique_filtered <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-Ub/processed_ratio_intensity/051313-tlj-52-55-matrix-ratios-ubi_unique_filtered.txt")

unfiltered = merge(vif_unfiltered, vpr_unfiltered, by=c("uniprot_id","Sequence"), all=T)
unfiltered = merge(unfiltered, vpu_unfiltered, by=c("uniprot_id","Sequence"), all=T)

ubi_filtered = merge(vif_ubifiltered, vpr_ubifiltered, by=c("uniprot_id","Sequence"), all=T)
ubi_filtered = merge(ubi_filtered, vpu_ubifiltered, by=c("uniprot_id","Sequence"), all=T)

ubi_unique_filtered = merge(vif_ubi_unique_filtered, vpr_ubi_unique_filtered, by=c("uniprot_id","Sequence"), all=T)
ubi_unique_filtered = merge(ubi_unique_filtered, vpu_ubi_unique_filtered, by=c("uniprot_id","Sequence"), all=T)

unfiltered_counts = apply(unfiltered[,3:8], 2, function(x) sum(!is.na(x)))
ubi_filtered_counts = apply(ubi_filtered[,3:8], 2, function(x) sum(!is.na(x)))
ubi_unique_filtered_counts = apply(ubi_unique_filtered[,3:8], 2, function(x) sum(!is.na(x)))

all_counts = cbind(unfiltered_counts, ubi_filtered_counts, ubi_unique_filtered_counts)
counts_selected_mean = t(apply(all_counts, 2,mean))
counts_selected_sd = t(apply(all_counts, 2, sd))
colnames(counts_selected_mean) = c("Total Peptides", "Ub Peptides", "Unique Peptides")

pdf("~/Projects/HPCKrogan/Data/HIV-proteomics/Plots/Ub-scans/Jurkat_filtering_steps.pdf", width=8, height=6)
barplot(t(all_counts), beside=T, col=c("lightgray","blue","gold"), ylab="Peptide counts", cex.names=0.8, legend.text=colnames(counts_selected_mean))
dev.off()

pdf("~/Projects/HPCKrogan/Data/HIV-proteomics/Plots/Ub-scans/Jurkat_filtering_steps_grouped.pdf", width=3, height=7)
barplot.EV(data=counts_selected_mean, errorvalues=counts_selected_sd,col=c("lightgray","blue","gold"), ylab="Peptide counts")
dev.off()
 