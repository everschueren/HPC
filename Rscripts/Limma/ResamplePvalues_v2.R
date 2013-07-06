library(sqldf)
library(stringr)
library(reshape)
source("Conversion/AnnotateWithUniprot.R")

## FUNCTIONS

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

## using ratios
input = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/293T-Expression-PTM/Vpu-Ub/processed/110112-Vpu-IFN-293T-Ub-matrix-max.txt", stringsAsFactors=F)

data_keys   = input[,1]
data_values = log2(input[,3:ncol(input)])
design = c(1,-1)
data_values = data_values * design
data_ranks  = convert_ratio_to_rank_matrix(data_values, normalized=F) 
data_ranked = cbind(data_keys, data_ranks) 
data_scored = score(data_ranked)
  
score = function(data_matrix){
  data_matrix = melt(data_matrix, id=c("data_keys"))[,c(1,3)]
  data_matrix = aggregate(. ~ data_keys, data=data_matrix, FUN=RP)  
}

sample_size=1000
PVALUE = 0.1

sample_scores = matrix(nrow=nrow(data_scored),ncol=sample_size, data=0)
for(i in 1:sample_size){
  ## reshuffle data_matrix
  for(c in 1:ncol(data_ranks)){
    data_ranks[,c] = sample(data_ranks[,c], size=nrow(data_ranks), replace=FALSE)
  }
  ## compute score
  tmp = cbind(data_keys, data_ranks)
  tmp_scored = score(tmp)
  ## attach score
  sample_scores[,i] = tmp_scored$value
}
## estimate p-value
p_values=c()
for(r in 1:nrow(sample_scores)){
  ecdf_sampled = ecdf(sample_scores[r,])
  key_score = ecdf_sampled(data_scored$value[r])
  p_values = c(p_values, key_score)
}
p_values = p.adjust(1-p_values, "BH")
data_scored$value = p_values
colnames(data_scored)[2] = "p_value"
data_scored_annotated = annotate_with_uniprot(data=data_scored, key="data_keys",organism="HUMAN")
data_scored_annotated_selected = sqldf(str_join("select * from data_scored_annotated where p_value < ", PVALUE, " order by p_value asc"))


## STOP HERE ###







#data_ranked_reshaped = sqldf("select * from data_ranked_reshaped order by data ASC")



## FLATTEN
# sample_names = colnames(data_ranks)
# sample_names = str_join(str_join("avg(",sample_names,")"),collapse=",")
# data_ranked_grouped = sqldf(str_join("select uniprot_id as 'uniprot_ac', ", sample_names,", count(*) as 'count' from data_ranked group by uniprot_id"))

PVALUE = 0.05
sample_size=1000
p_values = c()



## rank-resampling

data_ranks_avg = apply(data_ranks, 1, mean)
data_ranked = as.data.frame(cbind(data_keys, data_ranks_avg))
data_ranked = sqldf("select data_keys, avg(data_ranks_avg) as 'avg_rank' from data_ranked group by data_keys order by avg_rank asc")

for(i in 1:sample_size){
  # resample labels
  data_keys_reshuffled = sample(data_keys, size=length(data_keys), replace=FALSE)
  data_reshuffled = as.data.frame(cbind(data_keys_reshuffled, data_ranks_avg))
}





for(r in 1:nrow(input_values_grouped)){
  protein_id = input_values_grouped$uniprot_ac[r]
  
  ## get the number of peptides found for this proteins
  protein_peptides = input_values_grouped$count[r]
  
  ## get the avg. scores for this complex
  protein_score = input_values_grouped$average_score[r]
  sample_distro = c()
  
  ## make up random proteins of this size and compute their average 
  for(s in 1:sample_size){
    random_scores = sample(input_values_grouped$average_score, protein_peptides)
    sample_distro = c(sample_distro,mean(random_scores))
  }
  
  ## make a proper distribution of the sampled values
  ecdf_sampled = ecdf(sample_distro)
  
  ## compute the complex-score of the observed complex in the random sample ditribution
  ## the ecdf(c) function gives the number of observations equal or lower to x
  
  ecdf_score = ecdf_sampled(protein_score)
  if(protein_score>=0){
    p = 1-ecdf_score
  }else{
    p = ecdf_score
  }

  p_values = c(p_values,p)
  p_values_adj = p.adjust(p_values, "BH")
}

input_values_grouped = cbind(input_values_grouped, p_values_adj)
input_values_grouped_selected = sqldf(str_join("select * from input_values_grouped where average_score >= 0 and p_values_adj <= ",PVALUE," order by p_values_adj asc"))
#hist(input_values_grouped$p_values_adj)

input_values_grouped_selected = annotate_with_uniprot(data=input_values_grouped_selected, key="uniprot_ac",organism="HUMAN")