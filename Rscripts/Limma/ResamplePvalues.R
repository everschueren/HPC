library(sqldf)
library(stringr)
library(compiler)
library(limma)
library(gplots)
library(calibrate)
library(RColorBrewer)
source("Conversion/AnnotateWithUniprot.R")
source("Plotting/Aux.R")

## FUNCTIONS ##

resetPar = function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

## takes a frame [id, score] 
## computing the P that protein X with Y peptides has score Z  
resample.Peptides = function(scored_data, id_idx=1, score_idx=2, count_idx=3, sample_size=1000, two_sided=F){
  p_values = c()
  for(r in 1:nrow(scored_data)){
    
    id = scored_data[r,id_idx]
    ## get the avg. scores for this complex
    score = scored_data[r,score_idx]
    ## get the number of peptides found for this proteins
    count = scored_data[r,count_idx]
    
    sample_distro = c()
    
    ## make up random proteins of this size and compute their average 
    for(s in 1:sample_size){
      random_scores = sample(scored_data[,score_idx], count, replace=F)
      sample_distro = c(sample_distro,mean(random_scores))
    }
    
    ## make a proper distribution of the sampled values
    ecdf_sampled = ecdf(sample_distro)
    
    ## compute the protein-score of the observed group of peptides in the random sample ditribution
    ## the ecdf(c) function gives the number of observations equal or lower to x
    
    ecdf_score = ecdf_sampled(score)
    if(two_sided==F){
      p = 1-ecdf_score
    }else if(score<0 & two_sided==T){
      p = ecdf_score
    }
    p_values = c(p_values,p)
  }
  p_values_adj = p.adjust(p_values, "BH")
  p_values_adj
}
resample.Peptides.C = cmpfun(resample.Peptides)

##  DEVELOP
resample.Specificity = function(){
  
}

values_to_ranks = function(vec, normalized=T){
  tmp = as.numeric(factor(vec))
  m = max(tmp)
  if(normalized){
    (m - tmp + 1) / m  
  }else{
    (m - tmp + 1) 
  }
}

convert_ratio_to_rank_matrix = function(ratio_matrix, normalized=T){
  for(c in 1:ncol(ratio_matrix)){
    col = ratio_matrix[,c]
    ratio_matrix[,c] = values_to_ranks(col, normalized)
  }
  ratio_matrix
}

## compute Rank-Products
RP = function(vec){
  prod(vec)/length(vec)
}

## Compute dissimilarity
DS =  function(vec){
  mv = mean(vec)
  sum((vec - mv)^2)
}

resample.Reproducibility = function(scored_data, id_idx=1, score_idx=2, sample_size=1000, func=RP){
  ## convert scores to ranks
  data_ranks  = convert_ratio_to_rank_matrix(scored_data, normalized=F)
  data_RP = apply(data_ranks, 1, func)
  sampled_data_ranks = matrix(nrow=nrow(data_ranks), ncol=sample_size)
  for(s in 1:sample_size){
    sample_data_instance = matrix(nrow=nrow(data_ranks), ncol=ncol(data_ranks))
    for(col in 1:ncol(data_ranks)){
      sampled_col = sample(data_ranks[,col], size=nrow(data_ranks), replace=FALSE)
      sample_data_instance[, col] = sampled_col
    }
    sample_data_RP = apply(sample_data_instance, 1, func)
    sampled_data_ranks[, s] = sample_data_RP
  }
  p_values = c()
  for(d in 1:length(data_RP)){
    r_ecdf = ecdf(sampled_data_ranks[d,])
    r_pval = r_ecdf(data_RP[d])
    p_values = c(p_values, r_pval)
  }
  p_values_adj = p.adjust(p_values, "BH")
#   result = as.data.frame(cbind(1-(data_RP/max(data_RP))))
  result = as.data.frame(cbind(1-(data_RP/max(data_RP)), p_values_adj))
  colnames(result)[1] = "reproducibility_score" 
  colnames(result)[2] = "reproducibility_p_value"
  result
}

average.abundances = function(data, biolrep){
  b_leading = -9999
  abundance_avgs = c()
  idx = 1
  for(b in biolrep){
    if(b != b_leading){ # new biological sample
      
      #process previous sample, if any
      if(idx > 1){
        tmp = apply(abundance_tmp, 1, mean)
        abundance_avgs = cbind(abundance_avgs, tmp)
      }
        
      b_leading = b
      abundance_tmp = data[,idx]
    }else{ # same biological sample
      abundance_tmp = cbind(abundance_tmp, data[,idx])
    }
    idx = idx+1
  }
  ## process last one
  tmp = apply(abundance_tmp, 1, mean)
  abundance_avgs = cbind(abundance_avgs, tmp)
  colnames(abundance_avgs) = str_join("abundance_",rep(1:ncol(abundance_avgs)))
  abundance_avgs
}

########################
## INPUT & configuration

## general config
PVALUE = 0.05
SAMPLE_SIZE=100
ID_idx = 1
VAL_idx = 3

## specific config
dir = "~/Projects/HPCKrogan/Data/HIV-proteomics/Combined/processed/"
input_id = "combined-293"
# input_id = "combined-jurkat"
# input_id = "combined-deltas"
# input_id = "combined-wt"

NORMALIZE=T
PREPROCESSED_ABUNDANCES=F
PREPROCESSED_REPRODUCIBILITY=F

## READ
input_file = str_join(dir, input_id, "-matrix-max.txt")
input = read.delim(input_file)
id_name = colnames(input)[ID_idx]
## replace 1-ratios by NA to not screw up stats
input[input==1] = NA

#######################
## NORMALIZE INPUT DATA

if(NORMALIZE){
  tmp_data = as.matrix(input[,VAL_idx:ncol(input)])
  boxplot(log2(tmp_data),las=2, main="results before normalization", varwidth=T)
  tmp_data = normalizeBetweenArrays(tmp_data, method="scale")
  tmp_data = log2(tmp_data)
  boxplot(tmp_data,las=2, main="results after normalization", varwidth=T)
  input[,VAL_idx:ncol(input)] = tmp_data
}else{
  input[,VAL_idx:ncol(input)] = log2(input[,VAL_idx:ncol(input)]) 
}

##########################
## MAIN CALLS TO FUNCTIONS

####################################
## CONTINUE DEBUGGING SCRIPT HERE ##
####################################









###########################################################
## compute abundance ratios per protein from peptide ratios  

V_abundances = matrix(nrow=nrow(unique(input[ID_idx])),ncol=ncol(input)-VAL_idx+1)
colnames(V_abundances) = colnames(input)[VAL_idx:ncol(input)]
P_abundances = matrix(nrow=nrow(unique(input[ID_idx])),ncol=ncol(input)-VAL_idx+1)
colnames(P_abundances) = colnames(input)[VAL_idx:ncol(input)]

# par(oma=c(6,6,6,6))
# par(mar=c(5,4,4,2))
# heatmap.2(as.matrix(V_abundances[,1:24]))

abundance_file = str_join(dir, input_id, "-abundances.txt")
abundance_p_file = str_join(dir, input_id, "-abundances-pvalues.txt")

## if pre-processed then we read from file, else we compute the ratio-averages and p-values
if(PREPROCESSED_ABUNDANCES){
  P_abundances = read.delim(abundance_p_file, stringsAsFactors=F)
  V_abundances = read.delim(abundance_file, stringsAsFactors=F)
}else{
  ## cycle through all ratios
  for(idx in 1:ncol(V_abundances)){
    value_name = colnames(input)[VAL_idx+idx-1]
    ## compute statistics per protein
    input_grouped = sqldf(str_join("select ",id_name,", avg(",value_name,") as 'average_score', stdev(",value_name,") as 'stdev_score', count(*) as 'count' from input group by ",id_name," order by ",id_name," asc"))
    ## fix counts for proteins with 0 peptides (0=>1)
    input_grouped$count[is.na(input_grouped$average_score)] = 0 
    hist(input_grouped$average_score)
    
    P_abundance = resample.Peptides.C(scored_data=input_grouped, id_idx=1, score_idx=2, count_idx=3, sample_size=SAMPLE_SIZE)
    P_abundances[,idx] = P_abundance
    V_abundances[,idx] = input_grouped$average_score
  }
  ## write out processed values
  write.table(V_abundances, file=abundance_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  write.table(P_abundances, file=abundance_p_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
}

##########################
## compute reproducibility 

replicates = 2
replicate_design =c(-1,1)

exps = (ncol(V_abundances)/replicates)
V_reproducibility = matrix(ncol=exps,nrow=nrow(V_abundances), data=0)
P_reproducibility = matrix(ncol=exps,nrow=nrow(P_abundances), data=0)

reproducibility_file = str_join(dir, input_id, "-reproducibilities.txt")
reproducibility_p_file = str_join(dir, input_id, "-reproducibilities-pvalues.txt")

if(PREPROCESSED_REPRODUCIBILITY){
  V_reproducibility = read.delim(reproducibility_file)
  P_reproducibility = read.delim(reproducibility_p_file)
}else{
  for(i in 1:(ncol(V_abundances)/replicates)){
    i_start = ((i-1)*replicates)+1
    i_end = (i*replicates)
    abundance_data = V_abundances[,i_start:i_end] * replicate_design
    col_name = str_join(abundance_data,collapse="|")
    tmp = resample.Reproducibility(abundance_data, sample_size=1000, func=DS)
    V_reproducibility[,i] = as.matrix(tmp[,1])
    P_reproducibility[,i] = as.matrix(tmp[,2])
  }
  colnames(V_reproducibility) = rep("reproducibility",ncol(V_reproducibility))
  colnames(P_reproducibility) = rep("p_reproducibility",ncol(P_reproducibility))
  write.table(V_reproducibility,file=reproducibility_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  write.table(P_reproducibility,file=reproducibility_p_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)  
}

###########
## Specific

## specific config

vif_293 = 1
vpu_293 = 9
vpu_293_ifn = 10
vpu_293_ifn_mg132 = 11
vpu_293_mg132 = 12

vif_jurkat = 1
vpr_jurkat = 2
vpu_jurkat = 3

exp_idx = vif_jurkat
replicate_design =c(1,-1)

## GET ABUNDANCE
## work with replicates
# abundance_data =  V_abundances[,(((exp_idx-1)*replicates)+1):(exp_idx*replicates)]
## work with averaged replicates

abundance_data_repl = V_abundances[,(((exp_idx-1)*replicates)+1):(exp_idx*replicates)] 
abundance_data = average.abundances(abundance_data_repl * replicate_design, abs(replicate_design))

## GET REPRODUCIBILITY
## work with p-values
# abundance_data =  P_abundances[,exp_idx]
## work with scores
reproducibility_data = V_reproducibility[,exp_idx]
## work with p-values
# reproducibility_data = P_reproducibility[,exp_idx]

data_combined = cbind(input_grouped$uniprot_id, as.data.frame(abundance_data), as.data.frame(reproducibility_data))
colnames(data_combined) = c("uniprot_id", "abundance", "reproducibility")
data_combined_annotated = annotate_with_uniprot(data=data_combined, key="uniprot_id", organism="HUMAN")
data_combined_selected = sqldf("select * from data_combined_annotated where abundance != 0 order by abundance DESC")
plot(data_combined$abundance, data_combined$reproducibility, pch=20)
data_combined_selected[data_combined_selected$uniprot_id=="SCAM3_HUMAN",]
data_combined_selected[data_combined_selected$uniprot_id=="ITCH_HUMAN",]

