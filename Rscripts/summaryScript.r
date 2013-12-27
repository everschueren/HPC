#! /usr/bin/Rscript
suppressMessages(library(optparse))

# get the summarie values of each variable
getSums <- function(y){
	temp <- c(bait_name=unique(y$bait_name), bait_uniprot_id=unique(y$bait_uniprot_id), ID=unique(y$ID), num_prey=length(unique(y[,4])))
	temp <- c(temp, num_ranks=length(y$ms_rank[-grep("-", y$ms_rank)]))
	idx <- which(y[,6]==temp[2]) #location of bait in the prey list
	if(length(idx)>0){
		temp <- c(temp, bait_rank=y[idx,4], num_unique_peptides=y[idx,7], cov=y[idx,8])
	}else{
		temp <- c(temp, bait_rank="0", num_unique_peptides=0, cov=0)
	}
	temp <- c(temp, top_ranked_protein=y[which(y$ms_rank=="[1]"),]$ms_uniprot_ac)
	
	return(temp)
}


# finds which proteins are unique for that experiment
findProteins <- function(y, temp){
	prots <- which(table(y$ms_uniprot_ac)<2)
	idx <- match(names(prots), y$ms_uniprot_ac)
	tmp <- cbind(names(prots), y$ID[idx])
	tmp <- table(tmp[,2])
	tmp <- as.data.frame(cbind(names(tmp),tmp), stringsAsFactors=FALSE)
	names(tmp) <- c("ID", "original_peptides")
	tmp <- merge(temp, tmp, by="ID", all.x=TRUE)
	return(tmp)
}


# wrapper to summarize each bait ID
summarize <- function(x){
	tmp <- c()
	for(i in unique(x$ID)){
		idx <- which(x$ID == i)
		tmp <- rbind(tmp, getSums(x[idx,]))
	}
	tmp <- findProteins(x, tmp)
	# remove NA's
	tmp$original_peptides <- as.character(tmp$original_peptides)
	tmp$original_peptides[is.na(tmp$original_peptides)] <- 0
	tmp <- as.data.frame(tmp, stringsAsFactors=FALSE)
	tmp <- tmp[order(tmp$bait_rank),]
	return(tmp)
}

# get bait names and merge them with the Prospector data file
mergeWkeys <- function(datafile, keysfile){
	dat <- read.delim(datafile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
	keys <- read.delim(keysfile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
	names(dat)[1] <- names(keys)[2] <- "ID"
	dat <- merge(dat, keys[,c(2,5,7)], by="ID", all.x=TRUE)
	n <- dim(dat)[2]
	x1 <- dat[,c(n,n-1)]
	x2 <- dat[,-c(n,n-1)]
	dat <- cbind(x1, x2)
	return(dat)
}

# begin summarizing
main <- function(datafile, keysfile){
	dat <- mergeWkeys(datafile, keysfile)
	x <- summarize(dat)
	return(x)
}


option_list <- list(
	make_option(c("-d", "--data_file"), 
			help="data file containing Prospector output"),
	make_option(c("-o", "--output_file"), 
			help="output file for summary matrix"),
	make_option(c("-k", "--keys_file"), 
			help="keys file for summary matrix") #should be from MT          
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
summarizedData <- main(parsedArgs$data_file, parsedArgs$keys_file)

write.table(summarizedData, parsedArgs$output_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


###################################
#datafile="~/HPC/MSPipeline/tests/293T_VP/data/input/293T_VP_data.txt"
#keysfile="~/Desktop/12262013/KroganMS_2013_12_26_3-42 PM_MT_keys_merged.tsv"
#x <- main(datafile, keysfile)


