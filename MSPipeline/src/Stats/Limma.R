#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(limma))
suppressMessages(library(stringr))
suppressMessages(library(sqldf))
suppressMessages(library(optparse))
suppressMessages(library(compiler))
  
options(warn=-1)

# prepare_data = function(data_frame, convert_to_logs=T, normalize="none"){
#   data_matrix = as.matrix(data_frame[,3:ncol(data_frame)])
  
#   print(str_join("NORMALIZATION: ",normalize))
  
#   if(normalize!="none"){
#     data_matrix = normalizeBetweenArrays(data_matrix, method=normalize)
#   }
#   if(convert_to_logs){
#     data_matrix[data_matrix==0]=1
#     data_matrix = log2(data_matrix)
#   }
#   as.data.frame(data_matrix)
# }

make_design_vector = function(exp_samples, biol_replicates){
  rep(1:exp_samples, each=biol_replicates)
}

make_biol_tech_vector = function(exp_samples, biol_replicates){
  rep(1:(exp_samples*biol_replicates))
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

make_contrast_vector_from_file = function(contrast_file, sample_names){
  cf = as.data.frame(read.delim(contrast_file, header=F))
  contrasts = makeContrasts(contrasts=cf[,1], levels=sample_names) 
  contrasts
}

do_limma = function(data_matrix, design_matrix, contrasts=NULL){
  
  lin_fit = lmFit(data_matrix,design_matrix)  
  cnames = c("ID")
  
  if(is.null(contrasts)){
    eb = eBayes(lin_fit)
    tmp = topTable(eb, adjust.method="BH",number=Inf)
    tmp = cbind(rownames(tmp), tmp[,c(1,4)])
    exp_name = colnames(design_matrix)[1]
    colnames(tmp) =  c("ID","FC","Padj")
    results = tmp
  }else{
    contrasts.fit = eBayes(contrasts.fit(lin_fit, contrasts))
    #   test_results = decideTests(contrasts.fit, method="global")
    #   test_results_summary = summary(test_results)
    for(i in 1:ncol(contrasts)){
      test_differentials = topTable(contrasts.fit, coef=i, adjust.method="BH",number=Inf)
      tmp1 = cbind(rownames(test_differentials), test_differentials[,c(1,4)])
      ctr = colnames(contrasts)[i]
      colnames(tmp1)[1] = c("ID")
      
      if(i==1){
        results = tmp1
      }else{
        results = merge(results, tmp1, by="ID")
      }
      contrast_name = colnames(contrasts)[i]
      cnames = c(cnames, str_join(contrast_name,"_LFC"), str_join(contrast_name,"_Padj"))
    }
  }
  colnames(results) = cnames
  results
}

Limma.main = function(data_file, design_file, output_file, contrast_file="none"){
  data = read.delim(data_file, stringsAsFactors=F)
  design_matrix = read.delim(design_file)
  sample_names = colnames(design_matrix)
  data_matrix = data[,3:ncol(data)] 
  
  if(contrast_file=="none"){
    differential_results = do_limma(data_matrix, design_matrix) 
  }else{
    if(contrast_file=="auto"){
      print("AUTO CONTRASTS")
      contrasts = make_contrast_vector(sample_names=sample_names)
      print(contrasts)  
    }else{
      print(str_join("reading contrasts from file: ", contrast_file))
      contrasts = make_contrast_vector_from_file(contrast_file, sample_names)
    }
    differential_results = do_limma(data_matrix, design_matrix, contrasts) 
  }
    
  keys = cbind(rownames(data), data[,1:2])
  colnames(keys)[1] = "ID"
  output_frame = merge(keys, differential_results, by ="ID")
  write.table(output_frame, file=output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
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
  make_option(c("-m", "--design_file"),
              help="File containing design matrix for the experiment"),
  make_option(c("-c", "--contrast_file"), default="none",
              help="File containing the contrasts to be made from the design matrix")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
Limma.main <- cmpfun(Limma.main)
Limma.main(data_file=parsedArgs$data_file, design_file=parsedArgs$design_file, output_file=parsedArgs$output_file, contrast_file=parsedArgs$contrast_file)  

# Limma.main(data_file="~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-all-Ub/processed_repeats/Mock_v_WT_all_evidence_FLT_MAT.txt", design_file="~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-all-Ub/input/Mock_v_WT_all_design.txt", output_file="~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-all-Ub/processed_repeats/Mock_v_WT_all_evidence_LIM.txt", contrast_file="~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-all-Ub/input/Mock_v_WT_all_contrasts.txt")  




