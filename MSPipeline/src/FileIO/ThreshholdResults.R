#! /usr/bin/Rscript
suppressMessages(library(optparse))

main <- function(dat_file, comppass, mist, species){
	#dat_file="~/HPC/MSPipeline/tests/HHV8_glaunsinger/plots/293T_AndyVP_ALLSCORES.txt"
	x <- read.delim(dat_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)	
	comppass <- quantile(x$COMPPASS_WD, comppass)
	idx <- which(x$MIST_hiv>=as.numeric(mist) | x$COMPPASS_WD>=as.numeric(comppass))
	x <- x[idx,]
	
	species <- unlist(strsplit(species,"-"))
	idx <- which(!is.na(match(x$species, species)))
	return(x[idx,])	
}




option_list <- list(
	make_option(c("-d", "--data_file"), 
			help="data file containing Prospector output"),
	make_option(c("-o", "--output_file"), 
			help="output file for summary matrix"),
	make_option(c("-c", "--comppass_cutoff"), default=1,
			help="quantile of Comppass score to be used as a threshold"),
	make_option(c("-m", "--mist_cutoff"), default=.7,
			help="quantile of MiST score to be used as a threshold"),
		make_option(c("-s", "--species"),
			help="species we are interested in (removes others)")         
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))

filteredData <- main(parsedArgs$data_file, parsedArgs$comppass_cutoff, parsedArgs$mist_cutoff, parsedArgs$species)
write.table(filteredData, parsedArgs$output_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)




