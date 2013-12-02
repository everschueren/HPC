
#reads in a results file from MSPipeline
getScores <- function(filename){
	dat <- read.delim(filename, sep='\t', header=TRUE, stringsAsFactors =FALSE)
	dat <- dat2 <- dat[,c(1,2,3,5,11,14)]
	return(dat)
}


#create positive set from data&keys file
getPositives <- function(pos_data_file){
	print("Loading positive set")
	hnet <- read.delim(pos_data_file, sep="\t", header=TRUE, stringsAsFactors =FALSE)
	hnet <- hnet[hnet$score>2,]	#assume anything scored >2 is true interaction
	hnet <- hnet[,c(4,6)]		#get uniprot_ac's of pairs
	names(hnet) = c("bait_uniprot", "ms_uniprot_ac")
	return(unique(hnet))	
}




#get random sample of negatives: random pairing of all Human Genes (1% are true interactions)
getNegativeSample <- function(pos, human_gene_list){
	print("Getting True Negative set")
	#construct negative list
	#pos = getPositives(pos_data_file)
	mnegs <- read.delim(human_gene_list, sep='\t', header=TRUE, stringsAsFactors =FALSE)
	idx <- which(mnegs$Status == "reviewed")
	mnegs <- mnegs[idx,1]
	n <- dim(pos)[1]*99		#number of F interactions to T interactions ~ 99-1 ?
	idx <- sample(1:length(mnegs), n*2, replace=TRUE)
	negs = as.data.frame(matrix(mnegs[idx], ncol=2), stringsAsFactors=FALSE)
	names(negs) = c("bait_uniprot", "ms_uniprot_ac")
	
	#account for positive BAIT/PREY & PREY/BAIT scenario
	pos2 <- as.data.frame(cbind(bait_uniprot  = pos$ms_uniprot_ac, ms_uniprot_ac = pos$bait_uniprot), stringsAsFactors=FALSE)
	pos <-unique(rbind(pos,pos2))
	
	#make sure no duplicates b/w pos/neg sets
	truths = rbind(pos, negs)
	idx = which(duplicated(truths))
	while(length(idx)>0){	
		print(paste("REPLACING ", length(idx), " DUPLICATES", sep=""))
		n=length(idx)
		#get new sample to replace duplicates
		idx2 <- sample(1:length(mnegs), n*2, replace=TRUE)
		
		NEWnegs = matrix(mnegs[idx2], ncol=2)
		truths[idx,] = NEWnegs
		#idx <- which(duplicated(truths))
		idx <- checkCompliments(truths)
		idx <- idx[idx>dim(pos)[1]]		#dont replace the positives in the list
	}
	#return the list of non-duplicated negatives
	idx <- which(match(pos$bait_uniprot, truths$bait_uniprot) & match(pos$ms_uniprot_ac, truths$ms_uniprot_ac))
	negs = truths[-idx,]
	return(negs)	
}


checkCompliments <- function(x){
	#returns the index of any bait/prey = prey/bait pairings
	names(x) = c("a","b")
	x2 = cbind(b = x$a, a=x$b)
	x2 <- rbind(x,x2)
	idx <- sort(which(duplicated(x2)))
	idx1 <- idx[idx<=dim(x)[1]]
	idx2 <- idx[idx>dim(x)[1]] - dim(x)[1]
	return(unique(c(idx1,idx2)))	
}




#label truths as either positive= 1 OR negative= -1
scoreTruths <- function(pos, neg){
	print("Assigning Truth values")
	pos = cbind(pos, label=1)
	neg = cbind(neg, label = -1)
	x = rbind(pos, neg)
	names(x) = c("BAIT", "PREY", "label")
	return(x)
}


#create an index of all truth pairs to use in comparisons later (faster then direct comparison of pairs)
createIndex2 <- function(scores, truths){
	print("Indexing Bait/Prey Pairs")
	idx <- unique(rbind(truths[,1:2], scores[,1:2]))
	idx <- cbind(idx, idx = 1:dim(idx)[1])
	scores <- merge(scores, idx, by=c("BAIT", "PREY"), all.x=TRUE)
	truths <- merge(truths, idx, by=c("BAIT", "PREY"), all.x=TRUE)
	scores <- cbind(idx=scores[,dim(scores)[2]], scores[,-dim(scores)[2]])
	truths <- cbind(idx= truths[,dim(truths)[2]], truths[,-dim(truths)[2]])	
	return(list(scores, truths))
}




#create an index of all truth pairs to use in comparisons later (faster then direct comparison of pairs)
createIndex <- function(scores, truths){
	print("Indexing Bait/Prey Pairs")
	#must account for Bait/Prey and Prey/Bait matches
	scores2 <- as.data.frame(cbind(BAIT=scores$PREY, PREY=scores$BAIT), charactersAsFactors=FALSE)
	
	idx <- unique(rbind(truths[,1:2], scores[,1:2], scores2))	#5826665
	idx <- cbind(idx, idx = 1:dim(idx)[1])
	scores <- merge(scores, idx, by=c("BAIT", "PREY"), all.x=TRUE)
	scores2 <- merge(scores2, idx, by=c("BAIT", "PREY"), all.x=TRUE)
	truths <- merge(truths, idx, by=c("BAIT", "PREY"), all.x=TRUE)

	scores <- cbind(idx=scores[,dim(scores)[2]], scores[,-dim(scores)[2]])
	scores <- scores[order(scores$BAIT,scores$PREY),]
	scores2 <- scores2[order(scores2$PREY,scores2$BAIT),]	
	scores = merge(scores2, scores, by.x=c("PREY", "BAIT"), by.y=c("BAIT", "PREY"))

	truths <- cbind(idx= truths[,dim(truths)[2]], truths[,-dim(truths)[2]])	
	return(list(scores, truths))
}





#Run the respective pipeline score through with different thresholds to calculate ROC curve
getROCscores <- function(le_truths, le_scores){
	print(paste("Calculating ", names(le_scores)[2], " ROC scores"), sep="")
	tmp = data.frame()
	#cycle through threshhold values to create ROC curve
	for(thresh in 100:1){
		print("#############")
		print(thresh)
		
		x <- le_scores
		quant <- quantile(x[,3], thresh/100)
		idx1 <- which(x[,3] >= quant)
		idx2 <- which(x[,3] < quant)
		x[idx1,3] <- 1
		x[idx2,3] <- -1
		names(x)[3] <- "label"





		x = x[order(x$idx.x, x$idx.y),]
		head(x)
		x <- merge(x, truths[,c(1,4)], by.x="idx.y", by.y="idx", all.x=TRUE)	
		x = x[order(x$idx.x, x$idx.y),]
		head(x)
		x <- merge(x, truths[,c(1,4)], by.x="idx.x", by.y="idx", all.x=TRUE) #accounting for bait/prey or prey/bait possibility
		x = x[order(x$idx.x, x$idx.y),]
		head(x)
		
#!!!!!!!!!!!!!!!!!!!!!!!! checking on why label.y == label sometimes :(
		idx <- which(x$label == x$label.y)
		
		
		
		
		#fill in any values that were missed due to bait/prey or prey/bait ordering	
		idx = which(is.na(x$label.y) & !is.na(x$label))
		x$label.y[idx] <- x$label[idx]
		
		
		
		

		#x[is.na(x[,3]),] = -1	#anything that doesn't have a match = -1
		x = table(x[,3], x[,4])
		#if the predicted values predict all False or all True
		if(sum(dim(x) == c(2,2)) <2){
			if(rownames(x) == "-1"){
				tpr = fpr = 0	
			} else{
				tpr = fpr = 1
			}
		} else{
			x = x[2:1, 2:1]	#rotate confusion matrix
			print(x)
			tpr = x[1,1]/sum(x[,1])
			fpr = x[1,2]/sum(x[,2])
		}
		tmp = rbind(tmp, c(tpr, fpr))
		print(paste("tpr: ",round(tpr,4),"   fpr: ",round(fpr,4),sep=""))
	}
	names(tmp) = c("tpr", "fpr")
	return(tmp)
}


#generate ROC curves based on ROC scores
plotROC <- function(rocscores){
	plot(0,0, xlim=c(0,1), ylim=c(0,1), type='n', ylab="True Positive Rate", xlab="False Positive Rate", main="ROC Curve")
	for(i in 1:length(rocscores)){
		lines(rocscores[[i]]$tpr ~ rocscores[[i]]$fpr, col=i)
	}
	legend(.6,.2, names(rocscores), col=1:length(rocscores), lty=c(1,1), lwd=2.5)
}


#Main function that organizes scoring
scoreExperiment <- function(truths, scores){
	rocScores = list()
	#cycle through MIST_hiv, MIST_self, SAINT, etc scores
	for( scorecol in 5:dim(scores)[2] ){
		#go skim through multiple threshholds to create rock scores
		rocScores[[names(scores)[scorecol]]] = getROCscores(truths[,c(1,4)], scores[,c(3,4,scorecol)])
	}
	plotROC(rocScores)
	return(rocScores)
}







################
#-----MAIN------
################


#file that we will be using as benchmark
resultsfile = '~/HPC/MSPipeline/tests/humanNet/processed/humanNet_data_wKEYS_NoC_MAT_ALLSCORES.txt'
pos_data_file = '~/HPC/Benchmarks/datasets/humanNet/humannet_uniprot_keys.txt'
human_gene_list = '~/HPC/MSPipeline/files/uniprot_protein_descriptions_HUMAN.txt'	#used for negative list

scores = getScores(resultsfile)


#check to see if we alread have the Truths scored
truthsFile = "~/HPC/Benchmarks/datasets/humanNet/scoredTruths.txt"
if(file.exists(truthsFile)){
	print("Reading Pre-scored truths file")
	truths = read.delim(truthsFile, sep='\t', header=TRUE, stringsAsFactors =FALSE)
} else{
	print("Scoring truths and saving to file")
	pos = getPositives(pos_data_file)
	neg = getNegativeSample(pos, human_gene_list)
	#create list of truths based on our true positive and true negative sets
	truths <- scoreTruths(pos, neg)
	write.table(truths, truthsFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

#create an index of the pairs to hopefully speed up comparisons later
idx <- createIndex(scores, truths)
scores <- idx[[1]]
truths <- idx[[2]]

rocScores <- scoreExperiment(thruths, scores)














