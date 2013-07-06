library(stringr)
library(limma)
library(sqldf)

source("../Rscripts/Plotting/Aux.R")
source("../Rscripts/Plotting/MyPlots.R")

#######################################################
## FUNCTIONS ##########################################

## get counts without NA
getActualCounts = function(data){
  t(as.matrix(apply(data, 2, function(x) sum(!is.na(x)))))
}

## READ
readInput = function(input_file){
  input = read.delim(input_file)
  ## replace 1-ratios by NA to not screw up stats
  input[input==1] = NA
  input
}

## NORMALIZE 
normalize = function(data, method="scale"){ ## scale/cyclicloes/quantile
  tmp_data = as.matrix(data)
  tmp_data = normalizeBetweenArrays(tmp_data, method=NORMALIZE)
  tmp_data = log2(tmp_data)
  tmp_data
}

selectDFColsByRegEx = function(data, regex, invert=F, ignore.case=F){
  cnames = colnames(data)
  cols = grep(pattern=regex, x=cnames, value=F, invert=invert, ignore.case=ignore.case) 
  if(length(cols)==0){
    NULL
  }else{
    data[,cols] 
  }
}

na.dist <- function(x) {
  t.dist <- dist(x)
  t.dist <- as.matrix(t.dist)
  t.limit <- 1.1*max(t.dist,na.rm=T)
  t.dist[is.na(t.dist)] <- t.limit
  t.dist <- as.dist(t.dist)
  return(t.dist)
}

groups.cor = function(data){
  cor_data = cor(data, method="pearson", use="pairwise.complete.obs") 
  cor_data[1,2]
}

groups.labelflips = function(data){
  mask = apply(data, 1, function(x) sum(x) != sum(abs(x)))
  sum(mask, na.rm=T) / nrow(data)
}
# tmp = cbind(c(1,2,-1),c(1,-1,1))
# groups.labelflips(tmp)

groups.labels = function(data){
  str_join(colnames(data),collapse="\n")
}

groups.iterate = function(data, conditions, condition_names, proteins, protein_names, FUN, ACCFUN){
  accumulator = c()
  for(condition in conditions){
    condition_selected = selectDFColsByRegEx(data, condition, ignore.case=T)
    # print(colnames(condition_selected))
    if(is.null(condition_selected) == F){
      for(protein in proteins){
        protein_selected = selectDFColsByRegEx(condition_selected, protein, ignore.case=T)
        # print(colnames(protein_selected))
        if(is.null(protein_selected) == F){
          tmp = FUN(protein_selected)
          accumulator = ACCFUN(accumulator,tmp)
        }
      }  
    }
  }
  accumulator
}

#######################################################
## MAIN ###############################################

condition_names = c("-IFN/-MG132", "+IFN/-MG132", "-IFN/+MG132", "+IFN/+MG132")
conditions = c("[LH]$", "IFN$", "[LH]_MG132$", "IFN_MG132$")
proteins = c("GFP_[A-Z]_VIF|VIF_[A-Z]_GFP", "GFP_[A-Z]_VPR|VPR_[A-Z]_GFP", "GFP_[A-Z]_VPU|VPU_[A-Z]_GFP","GFP_[A-Z]_DeltaVIF|DeltaVIF_[A-Z]_GFP", "Mock_[A-Z]_WT|WT_[A-Z]_Mock")
protein_names = c("VIF", "VPR", "VPU","DeltaVIF", "MOCK")

ID_idx = 1
VAL_idx = 3
NORMALIZE = "scale"
# NORMALIZE = "cyclicloess"

dir = "~/Projects/HPCKrogan/Data/HIV-proteomics/Combined/processed/"
inputs = c("combined-293", "combined-jurkat", "combined-deltas", "combined-wt")
inputs = c("combined-293")

############# #####
## value collectors
counts = c()
col_names = c()
dup_cors_before = c()
dup_cors_after = c()
label_flips = c()
label_flips_afterN = c()
labels = c()

## MAIN DRIVER LOOP 

# TODO make plots decribing # of filtered out peptides 
# 1. (ub-sites vs non-ub)
# 2. (unique vs non-unique)

for(i in inputs){
  input_file = str_join(dir, i, "-matrix-max.txt")
  
  ## read input
  input = readInput(input_file)
  data = input[, VAL_idx:ncol(input)]
  keys = input[, ID_idx:(VAL_idx-1)]
  col_names = c(col_names, colnames(data))
  
  ## get counts
  tmp_counts = getActualCounts(data)
  counts = cbind(counts, tmp_counts)
  
  ## get silac-replication estimate
  # tmp = groups.iterate(data, conditions, condition_names, proteins, protein_names, groups.cor, c)
  # dup_cors_before = c(dup_cors_before, tmp)
  
  ## get selected labels
  tmp = groups.iterate(log2(data), conditions, condition_names, proteins, protein_names, groups.labels, c)
  labels = c(labels, tmp)
  
  ## get # of label-flips
  tmp = groups.iterate(log2(data), conditions, condition_names, proteins, protein_names, groups.labelflips, c)
  label_flips = c(label_flips, tmp)
  
  ## NEW PLOT: before boxplot
  # boxplot.EV(log2(data), main="not normalized", cex=0.6, xlab=input_id, ylab="Log2(SILAC-ratios)")
  
  ## normalize
  data_N = normalize(data, NORMALIZE)
  
  ## get silac-replication estimate and flipping label count after normalizations
  # tmp = groups.iterate(data_N, conditions, condition_names, proteins, protein_names, groups.cor, c)
  # dup_cors_after = c(dup_cors_after, tmp)
  
  ## NEW PLOT: after boxplot
  # boxplot.EV(data_N,main=str_join(NORMALIZE, " normalized"), cex=0.6, xlab=input_id, ylab="Log2(SILAC-ratios)")
  
  ## get # of label-flips
  tmp = groups.iterate(data_N, conditions, condition_names, proteins, protein_names, groups.labelflips, c)
  label_flips_afterN = c(label_flips_afterN, tmp)
  
  ## average per protein
  data_N_swapped = data_N * rep(c(1,1),ncol(data_N)/2)  ## re-swap label swaps (1,1) for identity or (-1,1) for label-reswapping
  data_N_swapped = as.data.frame(cbind(keys, data_N_swapped))
  data_N_grouped = sqldf(str_join("select uniprot_id, ", str_join(str_join("avg(",colnames(data),")"),collapse=","), " from data_N_swapped group by uniprot_id"))
  data_N_grouped = apply(as.matrix(data_N_grouped[,2:ncol(data_N_grouped)]),2,as.numeric)
  colnames(data_N_grouped) = colnames(data)

  ## NEW PLOT: protein heatmap
  # heatmap.2(t(data_N_grouped), distfun=na.dist, dendrogram="row",key=F,trace="none",col=bluered, margin=c(2,12))
  
  ## NEW PLOT: peptide heatmap
  # heatmap.2(t(data_N), distfun=na.dist, dendrogram="row",key=F,trace="none",col=bluered, margin=c(2,12))
  
  ## NEW PLOT: experiment heatmap
  ## use correlation for clustering
  cor_mat = cor(data_N_grouped, use="pairwise.complete.obs", method="pearson")
  # pdf(str_join("../../Data/HIV-proteomics/Plots/Ub-scans/normalized_experiment_clustering_",i,".pdf"))
  # heatmap.EV(cor_mat)
  # dev.off()
  
  ## NEW PLOT:
  ## protein overlap for conditions
  input_tmp =  input
  input_tmp[is.na(input_tmp)]=0 
  data_tmp = aggregate(. ~ uniprot_id, data=input_tmp,  FUN=mean)
  data_tmp[data_tmp==0]=NA
  overlaps = conditionOverlap(data_tmp[,VAL_idx:ncol(data_tmp)], conditions, condition_names, data_tmp[,1:(VAL_idx-1)])
  barplot(overlaps, beside=T, cex.names=0.6, ylim=c(0,4000),legend.text=c("condition 1", "condition 2", "overlap"))
}

mask.NArows = function(df){
  mask = apply(df,1,function(x)!all(is.na(x)))
  mask
}

# conditionOverlap = function(data, conditions, condition_names, keys){
#   data_tmp = data
#   data_tmp[is.na(data_tmp)] = 0
#   data_tmp[data_tmp!=0] = 1
#   data_tmp = as.matrix(data_tmp)
# }

conditionOverlap = function(data, conditions, condition_names, keys){
  res = c()
  cnames = c()
  for(i in 1:(length(conditions)-1)){
    condition1 = conditions[i]
    condition1_selected = selectDFColsByRegEx(data, condition1, ignore.case=T)
    c1_NNA_mask = mask.NArows(condition1_selected)
    c1_NNA = cbind(keys[c1_NNA_mask,], condition1_selected[c1_NNA_mask,])
    
    for(j in (i+1):length(conditions)){
      condition2 = conditions[j]
      if(condition1 != condition2){
        condition2_selected = selectDFColsByRegEx(data, condition2, ignore.case=T)
        c2_NNA_mask = mask.NArows(condition2_selected)
        c2_NNA = cbind(keys[c2_NNA_mask,], condition1_selected[c2_NNA_mask,])
        intersect = merge(c1_NNA, c2_NNA, by=c("uniprot_id"))
        print(str_join(c(condition1, condition2, nrow(c1_NNA), nrow(c2_NNA) ,nrow(intersect)), collapse=" "))
        res = cbind(res, c( nrow(c1_NNA), nrow(c2_NNA) ,nrow(intersect))) 
        cnames = c(cnames, str_join(condition_names[i],"\n",condition_names[j]))
      }
    }
  }
  colnames(res) = cnames
  res
}


##########
## NEW PLOT
## plot all Ub counts as barplot
counts = as.data.frame(counts)
barplot.EV(counts, cex=0.7, ylab="Ub sites", ylim_max=8000)

###########
## NEW PLOT
## pairplots estimating correlation between replicates and related experiments in single cell lines and condtions

label_flips = t(as.matrix(label_flips))
label_flips_afterN = t(as.matrix(label_flips_afterN))
colnames(label_flips) = labels
colnames(label_flips_afterN) = labels
barplot(label_flips, las=2,cex.names=0.5, ylim=c(0,1), main ="Pct. Before Normalization")
barplot(label_flips_afterN, las=2,cex.names=0.5, ylim=c(0,1), main ="Pct. After Normalization")

#####################################
## REWRITE CODE UNDER HERE AS GENERIC 

#########
## PLOT 5
## plot all Ub counts grouped per condition and cell type as barplot

# data_293 = data[,1:24]
# data_jurkat = data[,25:30]
# data_deltas = data[,31:34]
# data_wt = dat[,35:38]

##  estimate effect of IFN/MG132/BOTH on Ub counts
## group per MG132, IFN and w/o
# 
# conditiongroups = function(counts_current, conditions, condition_names){
#   means = c()
#   sds = c()
#   
#   for(condition in conditions){
#     counts_selected = selectDFColsByRegEx(counts_current, condition)
#     counts_selected_mean = apply(counts_selected,1,mean)  
#     counts_selected_sd = apply(counts_selected, 1, sd)
#     means = cbind(means, counts_selected_mean)
#     sds   = cbind(sds, counts_selected_sd)
#   }
#   
#   colnames(means) = condition_names
#   colnames(sds) = condition_names
#   list(means, sds)
# }
# 
# counts_293 = counts[,1:24]
# counts_jurkat = counts[,25:30]
# counts_deltas = counts[,31:34]
# counts_wt = counts[,35:38]
# 
# condition_names_293 = c("-IFN/-MG132", "+IFN/-MG132", "-IFN/+MG132", "+IFN/+MG132")
# conditions_293 = c("[LH]$", "IFN$", "[LH]_MG132$", "IFN_MG132$")
# condition_names_jurkat = c("-IFN/-MG132")
# conditions_jurkat = c("[LH]$")
# condition_names_deltas = c("-IFN/-MG132", "-IFN/+MG132")
# conditions_deltas = c("[LH]$","[LH]_MG132$")
# condition_names_wt = c("-IFN/-MG132","-IFN/+MG132")
# conditions_wt = c("[LH]$","[LH]_MG132$")
# proteins = c("VIF", "VPR", "VPU")
# protein_names = c("VIF", "VPR", "VPU")
# 
# mean_sds_293 = conditiongroups(counts_293, conditions_293, condition_names_293)
# # barplot.EV(mean_sds_293[[1]],errorvalues=mean_sds_293[[2]], ylim_max=8000, main="Ub sites 293 experiments")
# 
# mean_sds_jurkat = conditiongroups(counts_jurkat, conditions_jurkat, condition_names_jurkat)
# # barplot.EV(mean_sds_jurkat[[1]],errorvalues=mean_sds_jurkat[[2]], ylim_max=8000, main="Ub sites Jurkat experiments")
# 
# mean_sds_deltas = conditiongroups(counts_deltas, conditions_deltas, condition_names_deltas)
# # barplot.EV(mean_sds_deltas[[1]],errorvalues=mean_sds_deltas[[2]], ylim_max=8000, main="Ub sites Jurkat-∆ experiments")
# 
# mean_sds_wt = conditiongroups(counts_wt, conditions_wt, condition_names_wt)
# # barplot.EV(mean_sds_wt[[1]],errorvalues=mean_sds_wt[[2]], ylim_max=8000, main="Ub sites WT experiments")
# 
# ## all together
# barplot.EV(cbind(mean_sds_293[[1]],mean_sds_jurkat[[1]],mean_sds_deltas[[1]],mean_sds_wt[[1]]), errorvalues=cbind(mean_sds_293[[2]],mean_sds_jurkat[[2]],mean_sds_deltas[[2]],mean_sds_wt[[2]]),ylim_max=10000, main="Ub sites all experiments")
# 

