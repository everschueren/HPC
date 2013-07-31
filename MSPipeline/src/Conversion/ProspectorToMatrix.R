#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(reshape2))
suppressMessages(library(optparse))
suppressMessages(library(compiler))

options(warn=-1)

ProspectorToMatrix.removeIps = function(data_file, key_file){
  
} 

ProspectorToMatrix.cleanMatrix = function(data){
  ## get rid of specificity exclusion row and meta-columns
  data_tmp = data[-1,c(-2:-4)]
  colnames(data_tmp)[1] = "Preys"
  data_tmp
}

ProspectorToMatrix.cnames = function(data){
  cnames = colnames(data)
  cnames = cnames[c(-1:-4)]
  cnames
}

## dirty function to clean up matrix and convert it to long format
ProspectorToMatrix.convertMatrix = function(data, cnames){
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

ProspectorToMatrix.StatsTable = function(data_long){
  stats_tab = dcast(Preys ~ variable, data=data_long, fun.aggregate=mean, fill=0)
  rownames(stats_tab) = stats_tab[,1]
  stats_tab = as.matrix(stats_tab[,2:ncol(stats_tab)])
  stats_tab
}

ProspectorToMatrix.SpeciTable = function(stats_tab){
  stats_msk = stats_tab
  stats_msk[stats_msk>0]=1
  stats_msk
}

ProspectorToMatrix.ReproTable = function(data_long){
  repro_tab = dcast(Preys ~ variable, data=data_long, fun.aggregate=length, fill=0)
  rownames(repro_tab) = repro_tab[,1]
  repro_tab = as.matrix(repro_tab[,2:ncol(repro_tab)])
  repro_tab
}

ProspectorToMatrix.ResampledScreen = function(data_tmp){
  ## count # identified proteins in each run
  protein_cnts = apply(data.matrix(data_tmp[,2:ncol(data_tmp)]),2,function(x)sum(x>0))
  ## count # identified TSC in each run
  tsc_counts = apply(data.matrix(data_tmp[,2:ncol(data_tmp)]), 2, sum)
  ## cleanup
  Preys = data_tmp$Preys
  rownames(data_tmp) = Preys
  data_tmp = data_tmp[,-1]
  data_tmp = data.matrix(data_tmp, rownames.force=T)
  ## get TOTAL # of times prey identified
  prey_totals = apply(data_tmp, 1, sum)
  ## make random proteome 
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
    protein_total = protein_cnts[i]
    
    while(p < protein_total | tsc < tsc_total){
      
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
ProspectorToMatrix.ResampledScreen = cmpfun(ProspectorToMatrix.ResampledScreen, options=list("optimize",3))

data_file = "~/Projects/HPCKrogan/Data/HCV/Data/input/HCV-293T-All-results.txt"
master_file = "~/Projects/HPCKrogan/Data/Mastertable/processed/2013_May16_MasterTable_background_simple.txt"
remove_file = ""
keys_file = ""

ProspectorToMatrix.main = function(data_file, output_dir){
  print("READING")
  
  data = read.delim(data_file,stringsAsFactors=F, header=T,check.names=FALSE)
  
  ## remove ips
  
  ## merge with keys
  
  ## 
  
  if(master_file != ""){
    master =  read.delim(master_file, stringsAsFactors=F, header=T,check.names=FALSE)
  }
  
  
  print("CONVERTING")
  ## convert into intermediate format
  data_tmp = ProspectorToMatrix.cleanMatrix(data)
  cnames = ProspectorToMatrix.cnames(data)
  data_long = ProspectorToMatrix.convertMatrix(data_tmp, cnames)
  
  print("COMPUTING STATS")
  ## make stats table
  stats_tab = ProspectorToMatrix.StatsTable(data_long)
  
  ## make reproducibility table
  repro_tab = ProspectorToMatrix.ReproTable(data_long)
  
  ## make specificity table
  speci_tab = ProspectorToMatrix.SpeciTable(stats_tab)
    
  print("RESAMPLING SCREEN")
  ## resample the screen
  resampled_screen = ProspectorToMatrix.ResampledScreen(data_tmp)
  
  
}

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--data_file"),
              help="data file containing maxquant output"),
  make_option(c("-o", "--output_dir"),
              help="output file for converted matrix")
)

ProspectorToMatrix.main()

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
ProspectorToMatrix.main(parsedArgs$data_file,parsedArgs$output_dir)  


