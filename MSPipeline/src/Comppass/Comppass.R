library(reshape2)

data = read.delim("~/Projects/HPCKrogan/Scripts/MSPipeline/tests/comppass/HCV-293T-VP-results_NoC.txt", skip=1, stringsAsFactors=F)

data_good = grep(x=data$Rank, pattern='-',invert=T)
data = data[data_good,]

## STATS TABLE
stats_tab = dcast(Acc.. ~ X, data=data, value.var="Num.Unique", fun.aggregate= function(x) mean(x,na.rm=T), fill=0)
rownames(stats_tab) = stats_tab[,1]
stats_tab = as.matrix(stats_tab[,2:ncol(stats_tab)])

## STATS MASK [0=not found, 1=found]
stats_msk = stats_tab
stats_msk[stats_msk>0]=1

## REPRODUCIBILITY TABLE
repro_tab = dcast(Acc.. ~ X, data=data, value.var="Num.Unique", fun.aggregate=length, fill=0)
rownames(repro_tab) = repro_tab[,1]
repro_tab = as.matrix(repro_tab[,2:ncol(repro_tab)])


Comppass.Z = function(stats_tab){
  m = apply(stats_tab, 1, mean)
  sd = apply(stats_tab, 1, sd)
  Z = (stats_tab - m) / sd
  Z
}

Comppass.S = function(stats_tab, stats_msk){
  k = ncol(stats_msk)
  freq = apply(stats_msk, 1, sum)
  inv_freq = k / freq
  S = sqrt( inv_freq * stats_tab )
  S
}

Comppass.D = function(stats_tab, stats_msk, repro_tab){
  k = ncol(stats_msk)
  freq = apply(stats_msk, 1, sum)
  inv_freq = k / freq
  tmp = inv_freq * stats_tab
  
  for(i in 1:nrow(tmp)){
    for(j in 1:ncol(tmp)){
      tmp[i,j] = tmp[i,j] ^ (repro_tab[i,j])
    }
  }
  D = sqrt(tmp)
  D
}

write.table(D, file="~/Desktop/tmp.txt", sep="\t", quote=F, eol="\n", row.names=T, col.names=T)

Comppass.WD = function(){
  
}