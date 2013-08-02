#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(optparse))
suppressMessages(library(compiler))

INSTALL_DIR = Sys.getenv("MS_PIPELINE_PATH") 
source(paste(INSTALL_DIR,"src/Conversion/AnnotateWithUniprot_lib.R",sep=""))

CompileResults = function(dir="", output_file="", mist_metrics_file="", mist_self_score_file="", mist_hiv_score_file="", comppass_score_file="", saint_score_file="", filter_zeros=T, annotate=T, uniprot_dir="~/Projects/HPCKrogan/Scripts/MSPipeline/files/", species="HUMAN"){

	tmp = NULL
	header = c("BAIT", "PREY")

	if(mist_hiv_score_file!=""){
		mist_hiv_score = read.delim(paste(dir, mist_hiv_score_file, sep=""), stringsAsFactors=F)
		print("ADDING MIST_HIV_SCORES")
		tmp = mist_hiv_score
		header = c(header, c("MIST_hiv"))
	}

	if(mist_self_score_file!=""){
		mist_self_score = read.delim(paste(dir, mist_self_score_file, sep=""), stringsAsFactors=F)
		if(is.null(tmp)){
			tmp = mist_self_score
		}else{
			print("ADDING MIST_SELF_SCORES")
			tmp = merge(tmp, mist_self_score, by=c("Bait","Prey"))	
		}
		header = c(header, c("MIST_self"))
	}

	if(mist_metrics_file!=""){
		mist_metrics = read.delim(paste(dir, mist_metrics_file, sep=""), stringsAsFactors=F)
		if(is.null(tmp)){
			tmp = mist_metrics
		}else{
			print("ADDING MIST_METRICS")
			tmp = merge(tmp, mist_metrics, by=c("Bait","Prey"))	
		}
		header = c(header, c("MIST_R","MIST_A","MIST_S"))
	}

	if(comppass_score_file!=""){
		comppass_results = read.delim(paste(dir, comppass_score_file, sep=""), stringsAsFactors=F)
		if(is.null(tmp)){
			tmp = comppass_results
		}else{
			print("ADDING COMPPASS_SCORES")
			tmp = merge(tmp, comppass_results, by=c("Bait","Prey"))	
		}
		header = c(header, c("TSC_AVG","COMPPASS_Z","COMPPASS_S","COMPPASS_D","COMPPASS_WD","COMPPASS_pZ","COMPPASS_pS","COMPPASS_pD","COMPPASS_pWD"))
	}

	if(saint_score_file!=""){
		saint_results = read.delim(paste(dir, saint_score_file, sep=""), stringsAsFactors=F)
		if(is.null(tmp)){
			tmp = saint_results
		}else{
			print("ADDING SAINT_SCORES")
			tmp = merge(tmp, saint_results[,c("Bait","Prey","AvgP","MaxP")], by=c("Bait","Prey"))	
		}
		header = c(header, c("SAINT_AVG_P","SAINT_MAX_P"))
	}

	colnames(tmp) = header

	if(filter_zeros & comppass_score_file != ""){
		tmp = tmp[tmp$TSC_AVG > 0,]	
	}

	if(annotate){
		print("ANNOTATING WITH UNIPROT")
		# print(colnames(tmp))
		output = annotate_with_uniprot(as.data.frame(tmp), species=species, key="PREY", output_file=output_file, uniprot_dir=uniprot_dir)	
	}else{
		output = tmp
	}
	
}


option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--dir"),
              help="directory containing all the score files"),
  make_option(c("-o", "--output_file"),
              help="output file for converted matrix"),
  make_option(c("-t", "--mist_self_score_file"),
              help="file containing MIST self scores"),
  make_option(c("-f", "--mist_hiv_score_file"),
              help="file containing MIST hiv scores"),
  make_option(c("-m", "--mist_metrics_file"),
              help="file containing MIST metrics"),
  make_option(c("-c", "--comppass_score_file"),
              help="file containing COMPPASS scores"),
  make_option(c("-s", "--saint_score_file"),
              help="file containing SAINT scores"),
  make_option(c("-z", "--filter_zeros"), default=TRUE,
              help="filter unobserved bait-prey pairs from the scores"),
  make_option(c("-a", "--annotate"), default=TRUE,
              help="annotate results with uniprot fields"),
  make_option(c("-u", "--uniprot_dir"), 
              help="directory containing uniprot files"),
  make_option(c("-n", "--uniprot_species"), default="HUMAN", 
              help="species for which to read uniprot files"),
  make_option(c("-i", "--install_dir"), default="", 
              help="install dir of the pipeline")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))


CompileResults(parsedArgs$dir, parsedArgs$output_file, mist_metrics_file=parsedArgs$mist_metrics_file, mist_self_score_file=parsedArgs$mist_self_score_file, mist_hiv_score_file=parsedArgs$mist_hiv_score_file, comppass_score_file=parsedArgs$comppass_score_file, saint_score_file=parsedArgs$saint_score_file, filter_zeros=parsedArgs$filter_zeros, annotate=parsedArgs$annotate, uniprot_dir=parsedArgs$uniprot_dir, species=parsedArgs$uniprot_species)  


# CompileResults(dir="", output_file="/Users/everschueren/Projects/HPCKrogan/Data/HCV/Data/processed_v2/HCV-HuH-results_wKEYS_NoC_MAT_ALLSCORES.txt", mist_metrics_file="/Users/everschueren/Projects/HPCKrogan/Data/HCV/Data/processed_v2//MIST/HCV-HuH-results_wKEYS_NoC_MAT_MIST_SELF_metrics.txt", mist_self_score_file="/Users/everschueren/Projects/HPCKrogan/Data/HCV/Data/processed_v2//MIST/HCV-HuH-results_wKEYS_NoC_MAT_MIST_SELF_scores.txt", mist_hiv_score_file="/Users/everschueren/Projects/HPCKrogan/Data/HCV/Data/processed_v2//MIST/HCV-HuH-results_wKEYS_NoC_MAT_MIST_HIV_scores.txt", comppass_score_file="/Users/everschueren/Projects/HPCKrogan/Data/HCV/Data/processed_v2//COMPPASS/HCV-HuH-results_wKEYS_NoC_MAT_COMPPASS.txt", saint_score_file="/Users/everschueren/Projects/HPCKrogan/Data/HCV/Data/processed_v2/SAINT//RESULT/unique_interactions", filter_zeros=T, annotate=T, uniprot_dir="~/Projects/HPCKrogan/Scripts/MSPipeline/files/", species="HUMAN-HCV")
  