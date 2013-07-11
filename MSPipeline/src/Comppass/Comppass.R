#! /usr/bin/Rscript --vanilla --default-packages=utils
library(reshape2)
library(compiler)

Comppass.cleanMatrix = function(data){
  ## get rid of specificity exclusion row and meta-columns
  data_tmp = data[-1,c(-2:-4)]
  colnames(data_tmp)[1] = "Preys"
  data_tmp
}

Comppass.cnames = function(data){
  cnames = colnames(data)
  cnames = cnames[c(-1:-4)]
  cnames
}

## dirty function to clean up matrix and convert it to long format
Comppass.convertMatrix = function(data, cnames){
  ## make clean version of #replicated bait names w/o R re-naming
  cnames_rep = rep(cnames, each=nrow(data))
  ## convert into long format so we can later re-convert to wide format with an aggregate function
  data_l = melt(data, id=c("Preys"),  measure.vars=colnames(data)[2:ncol(data)])
  ## put the clean version of of #replicated bait names
  data_l$variable = cnames_rep
  ## make into numbers
  data_l$value = as.numeric(data_l$value)
  ## filter out 0's so that the mean gets computed correctly over # observations instead of # matrix cols   
  data_lf = data_l[data_l$value>0,]
  data_lf
}

Comppass.StatsTable = function(data_long){
  stats_tab = dcast(Preys ~ variable, data=data_long, fun.aggregate=mean, fill=0)
  rownames(stats_tab) = stats_tab[,1]
  stats_tab = as.matrix(stats_tab[,2:ncol(stats_tab)])
  stats_tab
}

Comppass.SpeciTable = function(stats_tab){
  stats_msk = stats_tab
  stats_msk[stats_msk>0]=1
  stats_msk
}

Comppass.ReproTable = function(data_long){
  repro_tab = dcast(Preys ~ variable, data=data_long, fun.aggregate=length, fill=0)
  rownames(repro_tab) = repro_tab[,1]
  repro_tab = as.matrix(repro_tab[,2:ncol(repro_tab)])
  repro_tab
}

Comppass.Z = function(stats_tab){
  m = apply(stats_tab, 1, mean)
  sd = apply(stats_tab, 1, sd)
  Z = (stats_tab - m) / sd
  Z
}

Comppass.S = function(stats_tab, speci_tab){
  k = ncol(speci_tab)
  freq = apply(speci_tab, 1, sum)
  inv_freq = k / freq
  S = sqrt( inv_freq * stats_tab )
  S
}

Comppass.D = function(stats_tab, speci_tab, repro_tab, normalized=F, D_T=1){ ## 95% score as normalization bar
  k = ncol(speci_tab)
  freq = apply(speci_tab, 1, sum)
  inv_freq = k / freq
  tmp = stats_tab
  
  for(i in 1:nrow(tmp)){
    for(j in 1:ncol(tmp)){
      tmp[i,j] = tmp[i,j] * (inv_freq[i] ^ (repro_tab[i,j]))
    }
  }
  D = sqrt(tmp)
  if(normalized){
    D_threshold = quantile(D, probs=D_T)
    D / D_threshold
  }else{
    D
  }
}

Comppass.WD = function(stats_tab, speci_tab, repro_tab, normalized=F, WD_T=1){ ## 95% score as normalization bar
  m = apply(stats_tab, 1, mean)
  sd = apply(stats_tab, 1, sd)
  w = sd / m
  w[w<=1]=1
  
  k = ncol(speci_tab)
  freq = apply(speci_tab, 1, sum)
  inv_freq = k / freq
  tmp = stats_tab
  
  for(i in 1:nrow(tmp)){
    for(j in 1:ncol(tmp)){
      tmp[i,j] = tmp[i,j] * (w[i] * inv_freq[i] ^ (repro_tab[i,j]))
    }
  }
  WD = sqrt(tmp)
  if(normalized){
    WD_threshold = quantile(WD, probs=WD_T)
    WD / WD_threshold
  }else{
    WD
  }
}

Comppass.Summary = function(stats_tab, Z, S, D, WD, WD_T=0){ ## WD_T takes everything over this  score
  
  score_table = c()
  
  for(i in 1:nrow(stats_tab)){
    for(j in 1:ncol(stats_tab)){
      if(WD[i,j] >= WD_T){
        bait = colnames(stats_tab)[j]
        prey = rownames(stats_tab)[i]
        z_score = Z[i,j]
        s_score = S[i,j]
        d_score = D[i,j]
        wd_score = WD[i,j]
        tsc = stats_tab[i,j]
        score_table = rbind(score_table, c(bait, prey, tsc, z_score, s_score, d_score, wd_score))  
      }
    }
  }
  score_table = as.data.frame(score_table)
  colnames(score_table) = c("Baits","Preys","Abundance","Z","S","D","WD")
  score_table
}
Comppass.Summary = cmpfun(Comppass.Summary, options=list("optimize",3))  

Comppass.ResampledScreen = function(data_tmp){
  protein_cnt = 300
  tsc_counts = apply(data.matrix(data_tmp[,2:ncol(data_tmp)]), 2, sum)
  Preys = data_tmp$Preys
  rownames(data_tmp) = Preys
  data_tmp = data_tmp[,-1]
  data_tmp = data.matrix(data_tmp, rownames.force=T)
  prey_totals = apply(data_tmp, 1, sum)
  random_proteome = rep(rownames(data_tmp), times=prey_totals)
  random_proteome_size = length(random_proteome)
  
  ## initialize random screen with 0 observations
  random_screen = data_tmp
  random_screen[] = 0
  
  for(i in 1:ncol(random_screen)){
    
    ## make a random run
    p=0
    tsc=0
    tsc_total = tsc_counts[i]
    
    while(p < protein_cnt & tsc < tsc_total){
      
      ## get a random protein and add it to run
      random_protein_idx = sample(1:random_proteome_size, 1)
      random_protein = random_proteome[random_protein_idx]
     
      current_tsc_count = random_screen[random_protein,i]
      if(current_tsc_count==0){ ## first time we sample this protein
        p = p + 1 
      }
      random_screen[random_protein,i] = current_tsc_count + 1
      tsc = tsc + 1
    }  
  }
  random_screen = cbind(Preys, random_screen)
  colnames(random_screen)[1] = "Preys"
  as.data.frame(random_screen)
}
Comppass.ResampledScreen = cmpfun(Comppass.ResampledScreen, options=list("optimize",3))


Comppass.main = function(data, output_file){
  print("CONVERTING")
  ## convert into intermediate format
  data_tmp = Comppass.cleanMatrix(data)
  cnames = Comppass.cnames(data)
  data_long = Comppass.convertMatrix(data_tmp, cnames)
  
  print("COMPUTING STATS")
  ## make stats table
  stats_tab = Comppass.StatsTable(data_long)
  
  ## make reproducibility table
  repro_tab = Comppass.ReproTable(data_long)
  
  ## make specificity table
  speci_tab = Comppass.SpeciTable(stats_tab)
  
  print("COMPUTING SCORES")
  ## compute scores
  Z = Comppass.Z(stats_tab)
  S = Comppass.S(stats_tab, speci_tab)
  D = Comppass.D(stats_tab, speci_tab, repro_tab)
  WD = Comppass.WD(stats_tab, speci_tab, repro_tab)
  WD_T = quantile(WD, probs=.95)
  
  print("RESAMPLING SCREEN")
  ## resample the screen
  random_screen = Comppass.ResampledScreen(data_tmp)
  random_data_long = Comppass.convertMatrix(random_screen, cnames)
  random_stats_tab = Comppass.StatsTable(random_data_long)
  random_repro_tab = Comppass.ReproTable(random_data_long)
  random_speci_tab = Comppass.SpeciTable(random_stats_tab)
  random_WD = Comppass.WD(random_stats_tab, random_repro_tab, random_speci_tab)
  random_WD_dist = ecdf(random_WD[random_WD>0])
  
  print("SUMMARIZING")
  ## compile all scores
  summary = Comppass.Summary(stats_tab, Z, S, D, WD, WD_T=0)
  summary = cbind(summary, 1.-random_WD_dist(summary$WD))
  colnames(summary)[8] = "P"
  
  print("WRITING")
  ## write out
  write.table(summary, file=output_file, row.names=F, col.names=T, eol="\n", sep="\t", quote=F) 
}

data = read.delim("~/Projects/HPCKrogan/Scripts/MSPipeline/tests/comppass/HCV-293T-Andy-results_wMT_MAT.txt",  stringsAsFactors=F, header=T,skip=1, check.names=FALSE)
Comppass.main(data, output_file="~/Projects/HPCKrogan/Scripts/MSPipeline/tests/comppass/HCV-293T-Andy-results_wMT_COMPPASS.txt")

data = read.delim("~/Projects/HPCKrogan/Scripts/MSPipeline/tests/comppass/HCV-293T-VP-results_wMT_MAT.txt",  stringsAsFactors=F, header=T,skip=1, check.names=FALSE)
Comppass.main(data, output_file="~/Projects/HPCKrogan/Scripts/MSPipeline/tests/comppass/HCV-293T-VP-results_wMT_COMPPASS.txt")

data = read.delim("~/Projects/HPCKrogan/Data/HCV/Data/processed/HCV-HuH-results_MAT.txt",  stringsAsFactors=F, header=T,skip=1, check.names=FALSE)
Comppass.main(data, output_file="~/Projects/HPCKrogan/Scripts/MSPipeline/tests/comppass/HCV-HuH-results_COMPPASS.txt")

# barplot(apply(data.matrix(data[2:nrow(data),5:ncol(data)]), 2, sum), las=2, cex.names=.3)
