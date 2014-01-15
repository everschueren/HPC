# This calculates the significance of overlap between two groups of proteins.
#	Using the classic urn approach, it assumes the whole human gene list as the population (20k+),
#	the first set A as the "white" set (Pass), the whole gene poulation minus this white set 
#	is our "black" set (Fail), and the second set B as the sampled set from
#	the whole gene list. 
# Specifically, these p-values are calculated on the overlap of each bait of groups A and B.
# This also applies a threshold based on MiST


####################
# calculates the set size, sample size, overlap, and significance (p-value) of overlap
overlapPerBait <- function(dat1, dat2, num_genes){
	# calculate overlap significance per bait of scored data
	temp <- c()
	baits <- unique(c(dat1$BAIT,dat2$BAIT))
	for(bait in baits){
		idx1 <- which(dat1$BAIT==bait)
		idx2 <- which(dat2$BAIT==bait)		
		# account for duplicate prey from consolidated datasets
		prey1 <- unique(dat1[idx1,]$PREY)
		prey2 <- unique(dat2[idx2,]$PREY)
		
		scores_o <- length(intersect(prey1, prey2))
		pval <- phyper(scores_o-1, length(prey1), num_genes-length(prey1), length(prey2), lower.tail=FALSE)
		tmp <- c(set=paste(unique(dat1$dataset), unique(dat2$dataset), sep="-"), BAIT=bait, num_preys1=length(prey1), num_preys2=length(prey2), num_overlap=scores_o, p_val=pval)
		temp <- rbind(temp, t(tmp))
		
	}
	
	# make the numbers numeric
	temp <- as.data.frame(temp, stringsAsFactors=FALSE)
	temp[,3] <- as.numeric(temp[,3])
	temp[,4] <- as.numeric(temp[,4])
	temp[,5] <- as.numeric(temp[,5])
	temp[,6] <- as.numeric(temp[,6])
	
	return(temp)
}

# inter-cell overlap
overlapPerCell <- function(dat1, dat2, num_genes, output_dir){
	d1 <- unique(dat1$PREY)
	d2 <- unique(dat2$PREY)
	
	#scores_o <- length(which(duplicated(c(d1, d2))))
	scores_o <- length(intersect(d1, d2))
	pval <- phyper(scores_o-1, length(d1), num_genes-length(d1), length(d2), lower.tail=FALSE)
	
	tmp <- c(set=paste(unique(dat1$dataset), unique(dat2$dataset), sep="-"), BAIT="_VIRUS_", num_preys1=length(d1), num_preys2=length(d2), num_overlap=scores_o, p_val=pval)
	
	interOverlapPreys(dat1, dat2, output_dir)
	return(tmp)
}

# spit out a file that labels overlapping and non-overlapping proteins
interOverlapPreys <- function(dat1, dat2, output_dir){
	d1 <- unique(dat1$PREY)
	d2 <- unique(dat2$PREY)
	overlaps <- intersect(d1, d2)
	set1 <- cbind(PREY=setdiff(d1, overlaps), dataset=unique(as.character(dat1$dataset)))
	set2 <- cbind(PREY=setdiff(d2, overlaps), dataset=unique(as.character(dat2$dataset)))
	dat_name <- paste(unique(set1[,2]), unique(set2[,2]), sep="-")
	overlaps <- rbind(set1, set2, cbind(PREY=overlaps, dataset = dat_name))[,2:1]
	
	if(!file.exists(paste( output_dir, "toEnrich", sep="")))
		dir.create(paste( output_dir, "toEnrich", sep=""))
	write.table(overlaps,paste(output_dir, "toEnrich/", dat_name, ".txt", sep=""),, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)	
}


#################
# wrapper to itterate through all combos of datasets and calculate the overlaps
getOverlap <- function(dat, thresh, output_dir){
	# normalize the BAIT/PREY names to all caps
	dat$BAIT <- toupper(dat$BAIT)
	dat$PREY <- toupper(dat$PREY)
	
	# find different datasets
	s <- as.character(unique(dat$dataset))
	comb <- t(combn(s,2))
	
	# set threshhold for data
	#thresh = 0.7
	dat <- dat[which(dat[,3] > thresh),]
	num_genes = 21121

	temp <- c()
	for( i in 1:dim(comb)[1]){
		dat1 <-  dat[which(dat$dataset == comb[i,1]),]
		dat2 <-  dat[which(dat$dataset == comb[i,2]),]
		if(unique(dat1$group)==unique(dat2$group)){
			x <- overlapPerBait(dat1, dat2, num_genes)	#calculates all the overlap info
		}else{
			x <- overlapPerCell(dat1, dat2, num_genes, output_dir)
		}
		temp <- rbind(temp, x)
	}
	return(temp)
}

# wrapper to just run and return the overlap and p-values
main <- function(dat, output_dir, thresh){
	if(!file.exists(output_dir))
		dir.create(output_dir)
		
	x <- getOverlap(dat, thresh, output_dir)
	write.table(x, paste(output_dir,"overlap.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	return(x)
}






###################
# read in data
#	This should be the only part that needs to be updated for future datasets
#	group = the case when we are comparing two separate viruses because they won't share baits
#-----------------------------
# read in the SCORES
#HHV8
andy <- cbind(read.delim("~/HPC/MSPipeline/tests/HHV8_glaunsinger/293T_Andy/data/processed/293T_Andy_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="293T", group=1)[,c(1,2,16,24,25)]
vp <- cbind(read.delim("~/HPC/MSPipeline/tests/HHV8_glaunsinger/293T_VP/data/processed/293T_VP_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="293T", group=1)[,c(1,2,16,24,25)]
bjab <- cbind(read.delim("~/HPC/MSPipeline/tests/HHV8_glaunsinger/BJAB/data/processed/BJAB_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="BJAB", group=1)[,c(1,2,16,24,25)]
islk <- cbind(read.delim("~/HPC/MSPipeline/tests/HHV8_glaunsinger/iSLK/data/processed/iSLK_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="iSLK", group=1)[,c(1,2,16,24,25)]

#HIV
hiv_hek <- cbind(read.delim("~/Box Documents/Projects/HIV/HIV_HEK.txt", sep="\t", stringsAsFactors=F), dataset="HIV_HEK", group=2)[,c(1,2,6,8,9)]
hiv_jurkat <- cbind(read.delim("~/Box Documents/Projects/HIV/HIV_Jurkat.txt", sep="\t", stringsAsFactors=F), dataset="HIV_JURKAT", group=2)[,c(1,2,6,8,9)]

#HCV
hcv_293andy <- cbind(read.delim("~/Box Documents/Projects/HCV/Data/processed_293Tandy_v7/293T_andy_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="HCV_293T", group=3)[,c(1,2,16,24,25)]
hcv_293vp <- cbind(read.delim("~/Box Documents/Projects/HCV/Data/processed_293Tvp_v7/293T_vp_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="HCV_293T", group=3)[,c(1,2,16,24,25)]
hcv_huh <- cbind(read.delim("~/Box Documents/Projects/HCV/Data/processed_huh_v7/huh_data_wKEYS_NoC_MAT_ALLSCORES.txt", sep="\t", stringsAsFactors=F), dataset="HCV_HUH", group=3)[,c(1,2,16,24,25)]

#scores <- rbind(andy, vp, islk, bjab, hiv_hek, hiv_jurkat, hcv_293andy, hcv_293vp, hcv_huh)
scores <- rbind(andy, vp, islk, bjab, hcv_293andy, hcv_293vp, hcv_huh)

##########################
output_dir <- "~/Desktop/overlaps/SAINT/"
#thresh = .7
thresh = quantile(scores[,3], 0.9)
x <- main(scores, output_dir, thresh)








