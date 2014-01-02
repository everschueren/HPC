#! /usr/bin/Rscript
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
##########################

makeQCPlots <- function(data_file, output_file){
	data_file_with_baits = data_file
	dat = read.delim(data_file_with_baits, header=TRUE, skip=1)
	names(dat) = c("X", "ID", "Rank", "Uniq.Pep", "Acc..", "Num.Unique", "X..Cov", "Best.Disc.Score", "Best.Expect.Val", "Protein.MW", "Species", "Protein.Name")
	bait_num = length(unique(dat[,1]))
	plots_per_col = 5
	plot_width = 15
	plot_height = ((plot_width/plots_per_col) * ceiling(bait_num/plots_per_col))
	# plot based on Num.Unique
	p = ggplot(dat, aes(x=factor(ID), y=Num.Unique)) ## consider using %coverage 
	
	pdf(file=paste(output_file, "Num_Unique_summary.pdf", sep="/"), width=plot_width, height=plot_height)
		print( p + geom_boxplot() + facet_wrap(~ X, scales="free", drop=T, ncol=plots_per_col) + theme(axis.text=element_text(size=9)) )
		## consider using scales=free_x to keep y scale fixed across all groups 
	dev.off()	
}
 

main <- function(data_file, output_file){
	makeQCPlots(data_file, output_file)
}

########################
option_list <- list(
	make_option(c("-d", "--data_file"),
              help="data file containing matrix"), #should be the "YOURDATANAME_data_wKEYS_NoC.txt" file in processed
	make_option(c("-o", "--output_file"),
              help="output directory for plots")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
main(parsedArgs$data_file, parsedArgs$output_file)




