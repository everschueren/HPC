#! /usr/bin/Rscript --vanilla --default-packages=utils

library(limma)
library(stringr)
library(sqldf)
library(optparse)
library(compiler)
scriptdir = "~/Projects/HPCKrogan/Scripts/Rscripts/"
source(str_join(scriptdir, "Conversion/AnnotateWithUniprot.R"))
  
options(warn=-1)

prepare_data = function(data_frame, convert_to_logs=T, normalize="none"){
  data_matrix = as.matrix(data_frame[,3:ncol(data_frame)])
  
  print(str_join("NORMALIZATION: ",normalize))
  
  if(normalize!="none"){
    data_matrix = normalizeBetweenArrays(data_matrix, method=normalize)
  }
  if(convert_to_logs){
    data_matrix[data_matrix==0]=1
    data_matrix = log2(data_matrix)
  }
  as.data.frame(data_matrix)
}

make_design_vector = function(exp_samples, biol_replicates, tech_replicates){
  rep(1:exp_samples, each=biol_replicates*tech_replicates)
}

make_biol_tech_vector = function(exp_samples, biol_replicates, tech_replicates){
  rep(1:(exp_samples*biol_replicates), each=tech_replicates)
}

make_sample_names = function(exp_samples){
  str_join("sample_", 1:exp_samples)
}

make_design_matrix = function(design_vector, sample_names){
  keys = unique(design_vector)
  dm = model.matrix(~0 + factor(design_vector))
  colnames(dm) = sample_names
  dm
}

make_contrast_vector = function(sample_names){
  contrasts = sample_names
  if(length(sample_names) >= 2){
    for(i in 1:(length(sample_names)-1)){
      for(j in (i+1):length(sample_names)){
        contrasts = c(contrasts, str_join(sample_names[i],"-",sample_names[j]))
      }
    } 
  }
  makeContrasts(contrasts=contrasts, levels=sample_names)
} 

make_contrast_vector_from_file = function(file, sample_names){
  contrast_file = read.delim(file)
  makeContrasts(contrasts=contrast_file$contrasts, levels=sample_names) 
}

do_limma = function(data_matrix, design_matrix, contrasts="none", biol_tech_vector=NULL, with_dup_cor=T){
  if(with_dup_cor ==T){
    cor_fit <- duplicateCorrelation(data_matrix,design_matrix,block=biol_tech_vector)
    print(str_join("cor.fit: ",cor_fit$cor))
  }
  if(with_dup_cor == T){ 
    if(is.nan(cor_fit$cor)== F & is.na(cor_fit$cor) == F){
      print("fitting with dup.cor")
      lin_fit = lmFit(data_matrix,design_matrix, block=biol_tech_vector, cor=cor_fit$cor)  
    }else{
      lin_fit = lmFit(data_matrix,design_matrix)  
    }
  }else{
    lin_fit = lmFit(data_matrix,design_matrix)  
  }
  if(contrasts == "none"){
    eb = eBayes(lin_fit)
    tmp = topTable(eb, adjust.method="BH",number=Inf)
    write.table(tmp, file="tmp.txt", eol="\n", sep="\t")
    tmp1 = cbind(rownames(tmp), tmp[,c(1,4)])
    exp_name = colnames(design_matrix)[1]
    colnames(tmp1) =  c("ID",str_join(exp_name,"_logFC"), str_join(exp_name,"_adjPVal"))
    tmp = tmp1
  }else{
    contrasts.fit = eBayes(contrasts.fit(lin_fit, contrasts))
    #   test_results = decideTests(contrasts.fit, method="global")
    #   test_results_summary = summary(test_results)
    for(i in 1:ncol(contrasts)){
      test_differentials = topTable(contrasts.fit, coef=i, adjust.method="BH",number=Inf)
      tmp1 = cbind(rownames(test_differentials), test_differentials[,c(1,5)])
      contrast = colnames(contrasts)[i]
      colnames(tmp1) = c("ID",str_join(contrast,"_logFC"), str_join(contrast,"_adjPVal"))
      
      if(i==1){
        tmp = tmp1
      }else{
        tmp = merge(tmp, tmp1, by="ID")
      }
    }
  }
  tmp
}

Limma.main = function(data_file, exp_samples, design_matrix_file="",normalize_method="none", contrast_file="none", biol_replicates, tech_replicates, output_file, convert_to_logs=T, with_dup_cor=T, organism="HUMAN"){
  data = read.delim(data_file, stringsAsFactors=F)
  data_matrix = prepare_data(data, convert_to_logs=T, normalize=normalize_method)
  biol_tech_vector = make_biol_tech_vector(exp_samples, biol_replicates, tech_replicates)

  if(design_matrix_file==""){
    sample_names = make_sample_names(exp_samples)
    design_vector = make_design_vector(exp_samples, biol_replicates, tech_replicates)
    design_matrix = make_design_matrix(design_vector, sample_names)  
  }else{
    print(str_join("reading design matrix from file: ", design_matrix_file))
    design_matrix = read.delim(design_matrix_file)
    sample_names = colnames(design_matrix)
  }
  if(contrast_file=="none"){
    contrasts = contrast_file
  }else if(contrast_file == "auto"){
    contrasts = make_contrast_vector(sample_names=sample_names)
    print("AUTO CONTRASTS")
    print(contrasts)
  }else{
    print(str_join("reading contrasts from file: ", contrast_file))
    contrasts = make_contrast_vector_from_file(contrast_file, sample_names)
  }
  
  differential_results = do_limma(data_matrix, design_matrix, contrasts, biol_tech_vector, with_dup_cor) 
  tmp = cbind(rownames(data), data[,1:2])
  colnames(tmp)[1] = "ID"
  output_frame = merge(tmp, differential_results, by ="ID")
  output_frame_annotated = annotate_with_uniprot(output_frame[,2:ncol(output_frame)], key="uniprot_id", organism=organism)
  write.table(output_frame_annotated,file=output_file,eol="\n",sep="\t",quote=F,row.names=F,col.names=T)
}


## READ OPTIONS

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-l", "--convert_to_logs"), default=TRUE,
              help="Flag to convert values in data-file to log2 values"),
  make_option(c("-d", "--data_file"),
              help="data file containing values"),
  make_option(c("-o", "--output_file"),
              help="output file for differential peptides"),
  make_option(c("-e", "--exp_samples"), type="integer",
              help="number of unique experimental sample sources"),
  make_option(c("-b", "--biol_reps"), type="integer", 
              help="number of biological replicates for sample sources"),
  make_option(c("-t", "--tech_reps"), type="integer", 
              help="number of technical replicates for sample sources"),
  make_option(c("-r", "--dup_cor"), default=TRUE,
              help="Flag to estimate duplicate correlation between probes"),
  make_option(c("-m", "--design_matrix"),
              help="File containing design matrix for the experiment"),
  make_option(c("-c", "--contrast_file"), default="none",
              help="File containing the contrasts to be made from the design matrix"),
  make_option(c("-n", "--normalize_method"), default="none",
              help="Normalization method: quantile, cyclicloess, none"),
  make_option(c("-u", "--uniprot_organism"), default="HUMAN",
              help="Uniprot species identifier to match hits to for more readable output")
)

## TEST ##############################################
TEST=FALSE

if(TEST){
  TESTDIR="~/Projects/HPCKrogan/Data/HIV-proteomics/Combined/"
  PREFIX="combined-293-VIF"
  NORMALIZATION="scale"
  REPLICACOMBINATION="max_ratio"
  TECH_REPLICAS=2
  BIOL_REPLICAS=2
  PVALUE=0.05
  MAXQCOLUMN="Ratio_H_L"
  FILTER="UBI"
  UNIQUE_FILTER=T
  DESIGN = "design_TBWT.txt"
  EXPERIMENTS=1
  ORGANISM="HUMAN"
  test_args = c("-d", str_join("processed/",PREFIX,"-matrix-max.txt"), "-o",str_join("processed/",PREFIX,"-differential.txt"), "-r", "FALSE", "-t", TECH_REPLICAS, "-e", EXPERIMENTS, "-b", BIOL_REPLICAS, "-m", str_join("processed/",DESIGN), "-n", NORMALIZATION, "-u", ORGANISM)
  setwd(TESTDIR)
  parsedArgs = parse_args(OptionParser(option_list = option_list), args = test_args)  
  debug(Limma.main)
  Limma.main(data_file=parsedArgs$data_file, design_matrix=parsedArgs$design_matrix, normalize_method=parsedArgs$normalize_method,contrast_file=parsedArgs$contrast_file,exp_samples=parsedArgs$exp_samples, biol_replicates=parsedArgs$biol_reps, tech_replicates=parsedArgs$tech_reps, output_file=parsedArgs$output_file, convert_to_logs=parsedArgs$convert_to_logs, with_dup_cor=parsedArgs$dup_cor, organism=parsedArgs$uniprot_organism)
}else{
  parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
  Limma.main_C <- cmpfun(Limma.main)
  Limma.main_C(data_file=parsedArgs$data_file, design_matrix=parsedArgs$design_matrix, normalize_method=parsedArgs$normalize_method,contrast_file=parsedArgs$contrast_file,exp_samples=parsedArgs$exp_samples, biol_replicates=parsedArgs$biol_reps, tech_replicates=parsedArgs$tech_reps, output_file=parsedArgs$output_file, convert_to_logs=parsedArgs$convert_to_logs, with_dup_cor=parsedArgs$dup_cor, organism=parsedArgs$uniprot_organism)  
}




