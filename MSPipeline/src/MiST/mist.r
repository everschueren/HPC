





processMatrix <- function(x){
	x <- x[-c(1:2),]
	names(x)[1:3] <- c("Preys", "PepAtlas", "Length")
	idx <- c(2,3,5:dim(x)[2])
	#convert columns from characters to numbers
	for(i in idx){
		x[,i] <- as.numeric(x[,i])
	}
	return(x)
}



getM3D <- function(x){
	x1 <- x[,5:dim(x)[2]]/x[,3]
	# normalize by column (of x)
	x1 <- scale(x1, center=FALSE, scale=colSums(x[,5:dim(x)[2]]))
	# normalize by column (of x1)	
	x1 <- scale(x1, center=FALSE, scale=colSums(x1))	
	# normalize by row
	temp <- t(scale(t(x1), center=FALSE, scale=colSums(t(x1))))

}














#dat <- read.csv("~/HPC/MSPipeline/tests/3A/input/3A-results.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#keys <- read.csv("~/HPC/MSPipeline/tests/3A/input/3A_keys.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

dat <- read.csv("~/HPC/MSPipeline/tests/3A/processed/3A-results_wKEYS_NoC_MAT.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

dat <- processMatrix(dat)
x <- dat

m3d <- getM3D(dat)








