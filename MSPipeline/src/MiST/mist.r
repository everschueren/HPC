




# simplify/convert data into workable form
processMatrix <- function(x){
	baits <- x[1,5:dim(x)[2]]
	x <- x[-c(1:2),]
	names(x)[1:3] <- c("Preys", "PepAtlas", "Length")
	idx <- c(2,3,5:dim(x)[2])
	#convert columns from characters to numbers
	for(i in idx){
		x[,i] <- as.numeric(x[,i])
	}
	return(list(x, baits))
}



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

}












##############################################################################################################
#dat <- read.csv("~/HPC/MSPipeline/tests/3A/input/3A-results.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#keys <- read.csv("~/HPC/MSPipeline/tests/3A/input/3A_keys.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

dat <- read.csv("~/HPC/MSPipeline/tests/3A_IP/processed/3A_IP-results_wKEYS_NoC_MAT.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

dat <- processMatrix(dat)
m3d_norm <- getM3D_normalized(dat[[1]])

#################
x <- m3d_norm
baits <- unlist(dat[[2]])

x <- getMetrics(m3d_norm, unlist(dat[[2]]))








