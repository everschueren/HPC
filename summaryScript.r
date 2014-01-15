#! /usr/bin/Rscript
suppressMessages(library(optparse))

# get the summarie values of each variable
getSums <- function(y){
	temp <- c(bait_name=unique(y$bait_name), bait_uniprot_id=unique(y$bait_uniprot_id), ID=unique(y$ID), num_unique_prey=length(unique(y[,4])))
	temp <- c(temp, num_unique_peptides=sum(na.omit(y$Uniq.Pep)))
	idx <- which(y$Acc..==temp[2]) #location of bait in the prey list
	if(length(idx)>0){
		temp <- c(temp, bait_detected="Y", bait_rank=y[idx,]$Rank)
	}else{
		temp <- c(temp, bait_detected="N", bait_rank="0")
	}
	temp <- c(temp, top_rank=y[which(y$Rank=="[1]"),]$Acc)
	
	return(temp)
}

# wrapper to summarize each bait ID
summarize <- function(x){
	tmp <- c()
	for(i in unique(x$ID)){
		idx <- which(x$ID == i)
		tmp <- rbind(tmp,getSums(x[idx,]))	
	}
	tmp <- as.data.frame(tmp)
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
			help="keys file for summary matrix")           
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
summarizedData <- main(parsedArgs$data_file, parsedArgs$keys_file)

write.table(summarizedData, parsedArgs$output_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)






