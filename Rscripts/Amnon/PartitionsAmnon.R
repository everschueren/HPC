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


## PARTITION SELECT FUNCTION
subsetSelectAmnon = function(machine, subset){
  res = NULL
  if(machine == "calcium"){
    if(subset == "infected"){
      ## Vif-interactors present in all Infected timepoints and absent in all Neg.Control timepoints (Calcium)
      res = prospectorSelectedPresence[prospectorSelectedPresence["AG_calcium-1"]==1 & prospectorSelectedPresence["AG_calcium-2"]==0 & prospectorSelectedPresence["AG_calcium-4"]==1 & prospectorSelectedPresence["AG_calcium-5"]==0 & prospectorSelectedPresence["AG_calcium-7"]==1 & prospectorSelectedPresence["AG_calcium-8"]==0 & prospectorSelectedPresence["AG_calcium-10"]==1 & prospectorSelectedPresence["AG_calcium-11"]==0,]  
    }else if(subset == "infected-stable"){
      ## Vif-interactors present in all Infected timepoints and absent in all Neg.Control timepoints and present in the Stable Cell line (Calcium)
      res = prospectorSelectedPresence[prospectorSelectedPresence["AG_calcium-1"]==1 & prospectorSelectedPresence["AG_calcium-2"]==0 & prospectorSelectedPresence["AG_calcium-3"]==1 & prospectorSelectedPresence["AG_calcium-4"]==1 & prospectorSelectedPresence["AG_calcium-5"]== 0 & prospectorSelectedPresence["AG_calcium-6"]==1 & prospectorSelectedPresence["AG_calcium-7"]==1 & prospectorSelectedPresence["AG_calcium-8"]==0 & prospectorSelectedPresence["AG_calcium-9"]==1 & prospectorSelectedPresence["AG_calcium-10"]==1 & prospectorSelectedPresence["AG_calcium-11"]==0 & prospectorSelectedPresence["AG_calcium-12"]==1 ,]
    }else if(subset == "infected-notstable"){
      ## Vif-interactors present in all Infected timepoints and absent in all Neg.Control timepoints and absent in the Stable Cell line (Calcium)
      res = prospectorSelectedPresence[prospectorSelectedPresence["AG_calcium-1"]==1 & prospectorSelectedPresence["AG_calcium-2"]==0 & prospectorSelectedPresence["AG_calcium-3"]==0 & prospectorSelectedPresence["AG_calcium-4"]==1 & prospectorSelectedPresence["AG_calcium-5"]==0 & prospectorSelectedPresence["AG_calcium-6"]==0 & prospectorSelectedPresence["AG_calcium-7"]==1 & prospectorSelectedPresence["AG_calcium-8"]==0 & prospectorSelectedPresence["AG_calcium-9"]==0 & prospectorSelectedPresence["AG_calcium-10"]==1 & prospectorSelectedPresence["AG_calcium-11"]==0 & prospectorSelectedPresence["AG_calcium-12"]==0 ,]
    }
  }else if(machine == "elite"){
    if(subset == "infected"){
      ## Vif-interactors present in all Infected timepoints and absent in all Neg.Control timepoints (Elite)
      res = prospectorSelectedPresence[prospectorSelectedPresence["AG-1"]==1 & prospectorSelectedPresence["AG-2"]==0 & prospectorSelectedPresence["AG-4"]==1 & prospectorSelectedPresence["AG-5"]==0 & prospectorSelectedPresence["AG-7"]==1 & prospectorSelectedPresence["AG-8"]==0 & prospectorSelectedPresence["AG-10"]==1 & prospectorSelectedPresence["AG-11"]==0,]
    }else if(subset == "infected-stable"){
      ## Vif-interactors present in all Infected timepoints and absent in all Neg.Control timepoints and present in the Stable Cell line (Calcium)
      res = prospectorSelectedPresence[prospectorSelectedPresence["AG-1"]==1 & prospectorSelectedPresence["AG-2"]==0 & prospectorSelectedPresence["AG-3"]==1 & prospectorSelectedPresence["AG-4"]==1 & prospectorSelectedPresence["AG-5"]== 0 & prospectorSelectedPresence["AG-6"]==1 & prospectorSelectedPresence["AG-7"]==1 & prospectorSelectedPresence["AG-8"]==0 & prospectorSelectedPresence["AG-9"]==1 & prospectorSelectedPresence["AG-10"]==1 & prospectorSelectedPresence["AG-11"]==0 & prospectorSelectedPresence["AG-12"]==1 ,]
    }else if(subset == "infected-notstable"){
      ## Vif-interactors present in all Infected timepoints and absent in all Neg.Control timepoints and absent in the Stable Cell line (Calcium)
      res = prospectorSelectedPresence[prospectorSelectedPresence["AG-1"]==1 & prospectorSelectedPresence["AG-2"]==0 & prospectorSelectedPresence["AG-3"]==0 & prospectorSelectedPresence["AG-4"]==1 & prospectorSelectedPresence["AG-5"]==0 & prospectorSelectedPresence["AG-6"]==0 & prospectorSelectedPresence["AG-7"]==1 & prospectorSelectedPresence["AG-8"]==0 & prospectorSelectedPresence["AG-9"]==0 & prospectorSelectedPresence["AG-10"]==1 & prospectorSelectedPresence["AG-11"]==0 & prospectorSelectedPresence["AG-12"]==0 ,]
    }
  }
  res = cbind(res, rep(str_join(machine, "_", subset),nrow(res)))
  colnames(res)[ncol(res)] ="subset"
  res
}

machine = "calcium"
machine = "elite"
subset = "infected"
subset = "infected-stable"
subset = "infected-notstable"

partition = subsetSelectAmnon("calcium","infected")
partition = subsetSelectAmnon("elite","infected")

quantifiedPartitions = merge(partition[,c(1,27)],prospectorSelectedValues,by="uniprot_id")
dcast(quantifiedPartitions, uniprot_id ~ subset , value.var=c('unique_peptides_num'), margins=T, fill=0,fun=mean)