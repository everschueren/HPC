# This calculates the significance of overlap between two groups of proteins.
#	Using the classic urn approach, it assumes the whole human gene list as the population,
#	the first set A as the "white" set (Pass), the whole gene poulation minus this white set 
#	is our "black" set (Fail), and the second set B as the sampled set from
#	the whole gene list.



####################
# calculates the set size, sample size, overlap, and significance (p-value) of overlap
overlapPerBait <- function(dat1, dat2, num_genes){
	# calculate overlap significance per bait of scored data
	temp <- c()
	baits <- unique(c(dat1$BAIT,dat2$BAIT))
	for(bait in baits){
		idx1 <- which(dat1$BAIT==bait)
		idx2 <- which(dat2$BAIT==bait)
		
		scores_o <- length(which(duplicated(c(scores1[idx1,]$PREY, scores2[idx2,]$PREY))))
		pval <- phyper(scores_o, length(idx1), num_genes-length(idx1), length(idx2), lower.tail=FALSE)
		
		tmp <- c(set=paste(unique(dat1$dataset), unique(dat2$dataset), sep="-"), num_preys1=length(idx1), num_preys2=length(idx2), num_overlap=scores_o, p_val=pval )
		temp <- rbind(temp, t(tmp))
		
	}
	
	# make the numbers numeric
	temp <- as.data.frame(temp, stringsAsFactors=FALSE)
	temp[,2] <- as.numeric(temp[,2])
	temp[,3] <- as.numeric(temp[,3])
	temp[,4] <- as.numeric(temp[,4])
	temp[,5] <- as.numeric(temp[,5])
	
	return(temp)
}


#################
# wrapper to itterate through all combos of datasets and calculate the overlaps
getOverlap <- function(dat, scores){
	# find different datasets
	s <- as.character(unique(dat$dataset))
	comb <- t(combn(s,2))
	
	# set threshhold for data
	thresh = 0.7
	dat <- dat[which(dat$MIST_self > thresh),]
	num_genes = 21121

	temp <- c()
	for( i in 1:dim(comb)[1]){
		dat1 <-  dat[which(dat$dataset == comb[i,1]),]
		dat2 <-  dat[which(dat$dataset == comb[i,2]),]
		x <- overlapPerBait(dat1, dat2, num_genes)	#calculates all the overlap info
		temp <- rbind(temp, x)
	}
	return(temp)
}



# wrapper to just run and return the overlap and p-values
main <- function(dat, output_file){
	x <- getOverlap(dat)	
	write.table(x, output_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	return(x)
}






###################
# read in data
#	This should be the only part that needs to be updated for future datasets
#-----------------------------
# read in the SCORES
andy <- cbind(read.delim("~/HPC/MSPipeline/tests/293T_Andy/data/processed/293T_Andy_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="293T")
vp <- cbind(read.delim("~/HPC/MSPipeline/tests/293T_VP/data/processed/293T_VP_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="293T")
bjab <- cbind(read.delim("~/HPC/MSPipeline/tests/BJAB/data/processed/BJAB_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="BJAB")
islk <- cbind(read.delim("~/HPC/MSPipeline/tests/iSLK/data/processed/iSLK_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="iSLK")
scores <- rbind(andy, islk, bjab)


output_file <- "~/Desktop/overlap.txt"


x <- main(scores, output_file)













