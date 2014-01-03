

# simplify/convert data into workable form
processMatrix <- function(x){
	baits <- x[1,5:dim(x)[2]]
	x <- x[-c(1:2),]
	names(x)[1:3] <- c("Preys", "PepAtlas", "Length")
	row.names(x) <- x$Preys
	idx <- c(2,3,5:dim(x)[2])
	#convert columns from characters to numbers
	for(i in idx){
		x[,i] <- as.numeric(x[,i])
	}
	return(list(x, baits))
}

# "normalize" the M3D variable
getM3D_normalized <- function(x){
	x1 <- x[,5:dim(x)[2]]/x[,3]
	# normalize by column (of x)	-	M3D
	x1 <- scale(x1, center=FALSE, scale=colSums(x[,5:dim(x)[2]]))
	# normalize by column (of x1)	-	M3D normalized
	x1 <- scale(x1, center=FALSE, scale=colSums(x1))	
	return(x1)
}

# calculate Abundance, Reproducibility, and Specificity
getMetrics <- function(x, baits){
	reproducibility <- c()
	abundance <- c()
	
	# look at data per bait
	for( i in unique(baits)){
		bidx <- which(baits == i)
		y <- x[,bidx]
		# normalize by row by bait
		y <- y/apply(y, 1, sum)
		y[is.na(y)] <- 0
		
		# Reproducibility ("entropies")
		# -----------------------------
		y[y!=1 & y>0] = y[y!=1 & y>0] * log2(y[y!=1 & y>0])
		y[y==1] = (1-1e-10)*log2(1-1e-10)
		y <- rowSums(y)
		if(length(bidx)!=1){
			y = y/log2(1/length(bidx))
		}else{
			y = y*0+1
		}
		reproducibility <- cbind(reproducibility, y)
		colnames(reproducibility)[dim(reproducibility)[2]] <- i
		
		# Abundance ("averages")
		# ----------------------
		abundance <- cbind(abundance, rowSums(x[,bidx])/length(bidx))
		colnames(abundance)[dim(abundance)[2]] <- i
	}
	
	# Specificity
	# ----------------------
	specificity <- t(apply(abundance,1,function(z) z/sum(z)))
	
	return(list(reproducibility, abundance, specificity))
}

# vectorize the metrics while keeping the names straight
vectorize <- function(x){
	temp <- c()
	for(i in 1:dim(x)[2] ){
		tmp <- as.data.frame(cbind(x[,i], colnames(x)[i], names(x[,i])), stringsAsFactors=FALSE)
		temp <- rbind(temp,tmp)
	}
	names(temp) <- c("Xscore", "Bait", "Prey")
	temp$Xscore <- as.numeric(temp$Xscore)
	return(temp[,c(3,2,1)])
}


# Perform PCA analysis AS DONE IN MIST
doPCA <- function(R,A,S){
	# vectorize the variables -- may need to make changes to carry over names
	m <- cbind(R=R$Xscore, A=A$Xscore, S=S$Xscore )
	x <- princomp(m)
	
	#now do some other stuff??
	scores <- -x$scores[,1]
	scores <- 1 - (scores - min(scores))/(max(scores)-min(scores))
	scores <- cbind(R=R, A=A$Xscore, S=S$Xscore, MiST=scores)
	names(scores) = c("Prey", "Bait", "Reproducibility", "Abundance", "Specificity", "MiST_score")
	return(scores)
}


##############################################################################################################
#dat <- read.csv("~/HPC/MSPipeline/tests/3A/input/3A-results.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#keys <- read.csv("~/HPC/MSPipeline/tests/3A/input/3A_keys.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

#dat <- read.csv("~/HPC/MSPipeline/tests/3A_IP/processed/3A_IP-results_wKEYS_NoC_MAT.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
dat <- read.csv("~/HPC/MSPipeline/tests/HHV8_glaunsinger/BJAB/data/processed/BJAB_data_wKEYS_NoC_MAT.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)


dat <- processMatrix(dat)
m3d_norm <- getM3D_normalized(dat[[1]])

dat <- getMetrics(m3d_norm, unlist(dat[[2]]))
R <- vectorize(dat[[1]])
A <- vectorize(dat[[2]])
S <- vectorize(dat[[3]])

x <- doPCA(R,A,S)

write.table(x, "~/HPC/MSPipeline/tests/HHV8_glaunsinger/BJAB/data/processed/MiST_scores_1.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t" )




