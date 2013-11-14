
#reads in a results file from MSPipeline
getScores <- function(filename){
	dat <- read.delim(filename, sep='\t', header=TRUE, stringsAsFactors =FALSE)
	dat <- dat2 <- dat[,c(1,2,3,5,11,14)]
	#now account for bait/pret and prey/bait combos
	#names(dat2)[1:2] = c("PREY","BAIT")
	#dat <- rbind(dat, dat2)
	return(dat)
}


#create positive set from data&keys file
getPositives <- function(){
	print("Getting True Positve set")
	#construct positivelist
	pos <- read.delim("~/HPC/Benchmarks/datasets/humanNet/humanNet.txt", sep="\t", header=TRUE, stringsAsFactors =FALSE)
	keys <- read.delim('~/HPC/Benchmarks/datasets/humanNet/humanNet_keys.txt', sep='\t', header=TRUE, stringsAsFactors =FALSE)
	dat <- merge(pos[,c(1,4)], keys, by="id")[,c(3,2)]
	return(unique(dat))	
}


#get random sample of negatives: random pairing of all Human Genes (1% are true interactions)
getNegativeSample <- function(s, pos){
	print("Getting True Negative set")
	#construct negative list
	mnegs <- read.delim('~/HPC/MSPipeline/files/uniprot_protein_descriptions_HUMAN.txt', sep='\t', header=TRUE, stringsAsFactors =FALSE)
	idx <- which(mnegs$Status == "reviewed")
	mnegs <- mnegs[idx,1]	
	idx <- sample(1:length(mnegs), s*2, replace=TRUE)
	negs = as.data.frame(matrix(mnegs[idx], ncol=2))
	names(negs) = c("bait_uniprot", "ms_uniprot_ac")
	
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
		idx = which(duplicated(truths))
	}
	
	#return the list of non-duplicated negatives
	idx <- which( (pos[,1]==truths[,1]) & (pos[,2]==truths[,2]) )
	negs = truths[-idx,]
	return(negs)	
}


#label truths as either positive= 1 OR negative= -1
scoreTruths <- function(pos, neg){
	print("Assigning Truth values")
	pos = cbind(pos, label=1)
	neg = cbind(neg, label = -1)
	x = rbind(pos, neg)
	names(x) = c("BAIT", "PREY", "label")
	return(unique(x))
}


#create an index of all truth pairs to use in comparisons later (faster then direct comparison of pairs)
createIndex <- function(scores, truths){
	print("Indexing Bait/Prey Pairs")
	idx <- unique(rbind(truths[,1:2], scores[,1:2]))
	idx <- cbind(idx, idx = 1:dim(idx)[1])
	scores <- merge(scores, idx, by=c("BAIT", "PREY"), all.x=TRUE)
	truths <- merge(truths, idx, by=c("BAIT", "PREY"), all.x=TRUE)
	scores <- cbind(idx=scores[,dim(scores)[2]], scores[,-dim(scores)[2]])
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
		
		x = le_scores
		quant = quantile(x[,2], thresh/100)
		idx1 = which(x[,2] >= quant)
		idx2 = which(x[,2] < quant)
		x[idx1,2] = 1
		x[idx2,2] = -1
		names(x)[2] = "label"
		
		x = merge(x, truths[,c(1,4)], by = c("idx"), all.x=TRUE)
		#x[is.na(x[,3]),] = -1	#anything that doesn't have a match = -1
		x = table(x[,2], x[,3])
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
	for( scorecol in 4:dim(scores)[2] ){
		#go skim through multiple threshholds to create rock scores
		rocScores[[names(scores)[scorecol]]] = getROCscores(truths[,c(1,4)], scores[,c(1,scorecol)])
	}
	plotROC(rocScores)
	return(rocScores)
}







################
#-----MAIN------
################


#file that we will be using as benchmark
resultsfile = '~/HPC/MSPipeline/tests/humanNet/processed/humanNet_data_wKEYS_NoC_MAT_ALLSCORES.txt'
scores = getScores(resultsfile)
#scores <- cbind(scores, invMIST=(1-scores$MIST_self) )
#scores <- scores[,c(1,2,4)]

#check to see if we alread have the Truths scored
truthsFile = "~/HPC/Benchmarks/datasets/humanNet/scoredTruths.txt"
if(file.exists(truthsFile)){
	print("Reading Pre-scored truths file")
	truths = read.delim(truthsFile, sep='\t', header=TRUE, stringsAsFactors =FALSE)
} else{
	print("Scoring truths and savind to file")
	pos = getPositives()
	neg = getNegativeSample(dim(pos)[1]*99, pos)
	#create list of truths based on our true positive and true negative sets
	truths <- scoreTruths(pos, neg)
	write.table(truths, truthsFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

#create an index of the pairs to hopefully speed up comparisons later
idx <- createIndex(scores, truths)
scores <- idx[[1]]
truths <- idx[[2]]

rocScores <- scoreExperiment(thruths, scores)














