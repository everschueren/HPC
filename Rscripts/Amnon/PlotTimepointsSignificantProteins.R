library(gplots)
library(stringr)

print_heatmap = function(timepoints_file, output_file, main_suffix=""){
  timepoints = read.delim(timepoints_file)
  input_matrix = as.matrix(timepoints[,2:5])
  rownames(input_matrix) = timepoints[,1]  
  rc <- rainbow(nrow(input_matrix), start = .0, end = .3)
  pdf(width=10,height=10,file=output_file)
  bname = basename(timepoints_file)
  heatmap.2(input_matrix, Colv = NA, scale = "none", col=rc, margins = c(12,12), cexRow = 0.1 + 1/(log10(nrow(input_matrix))+2), dendrogram='row',keysize = 1.5,density.info="none", trace="none", main=str_join("Log2-intensities timepoints ",main_suffix))  
  dev.off()
}

print_heatmap("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PPI/processed/011013-ag-1-12-differential_infected_neg_control_selected_timepoints.txt", "~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PPI/processed/011013-ag-1-12-differential_infected_neg_control_selected_timepoints.pdf", main="Vif (round 0)")
print_heatmap("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PPI/processed/041213-ag-13-24-differential_infected_neg_control_selected_timepoints.txt", "~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PPI/processed/041213-ag-13-24-differential_infected_neg_control_selected_timepoints.pdf", main="Vif (round 1)")
print_heatmap("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PPI/processed/041313-ag-25-36-differential_infected_neg_control_selected_timepoints.txt", "~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PPI/processed/041313-ag-25-36-differential_infected_neg_control_selected_timepoints.pdf", main="Vpr (round 1)")
print_heatmap("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PPI/processed/041513-ag-37-48-differential_infected_neg_control_selected_timepoints.txt", "~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PPI/processed/041513-ag-37-48-differential_infected_neg_control_selected_timepoints.pdf", main="Vpu (round 1)")