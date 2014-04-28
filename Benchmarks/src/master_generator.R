

###############################
# Formatting HNet data
#	Find which hnet data we have in our MT then save data and keys file
dat <- read.delim("~/HPC/Benchmarks/datasets/master/MT_data_20131024.txt", sep="\t", header=TRUE, stringsAsFactors =FALSE)
keys <- read.delim("~/HPC/Benchmarks/datasets/master/MT_keys_20131024.txt", sep="\t", header=TRUE, stringsAsFactors =FALSE)

#create useable master bait/prey table
master = merge(keys[,c(1,9)], dat, by="id")
#Filter out un-usefull entries
master <- master[-which(master$ms_uniprot_ac == 'decoy'),]
master <-  master[-which(as.numeric(master$ms_uniq_pep) <= 0),]	#NA warning b/c "" 's are turned into NA's. Is ok


hnet <- read.delim('~/Dropbox/kroganlab/dataset/apms_benchmark/datasets/HNet/humannet_uniprot_keys.txt', sep='\t', header=TRUE, stringsAsFactors =FALSE)
#keep only scores > 2
hnet <- hnet[hnet$score>2,]


#create pairs to compare
pears = list()
pears[["pair1"]] = paste(hnet$uniprot_ac, hnet $uniprot_ac.1, sep='-')
pears[["pair2"]] = paste(hnet $uniprot_ac.1, hnet $uniprot_ac, sep='-')
pears[["m1"]] = paste(master$bait_uniprot, master$ms_uniprot_ac, sep="-")


#Function that returns the 
comparePairs <- function(pairSet1, pairSet2){
	#pull all the baits that are common to this pairing
	idx <- which(!is.na(match(pairSet2, pairSet1)))	#returns location of matched pairs from the master list
	print(length(idx))
	x = unique(unlist(lapply(strsplit(pairSet2[idx],"-"), function(x) x[1])))
	print(length(x))
	return(list(x, idx))
}

#compile list of baits to pull from the masterset
#	save their idx from master set in order to create keys file
toget <- comparePairs(pears[[1]], pears[[3]])
toget <- c(toget, comparePairs(pears[[2]], pears[[3]]) )

idx <- unique(c(toget[[2]], toget[[4]]))
toget = sort(unique(c(toget[[1]], toget[[3]])), decreasing=TRUE)

#pull all baits from the masterset
hnet <- master[which(!is.na(match(master$bait_uniprot, toget))),]

#write out the hnet positive set with threshold
write.table(as.matrix(hnet[,-2]), "~/HPC/Benchmarks/datasets/humanNet/humanNet_data.txt", sep="\t", row.names=FALSE, quote=FALSE)


######
#write the keys file for the hnet file
x = unique(hnet[,c(1,2)])
write.table(as.matrix(x), "~/HPC/Benchmarks/datasets/humanNet/humanNet_keys.txt", sep="\t", row.names=FALSE, quote=FALSE)


#clean up memory (remove variables)
rm(list=ls())




###############################
#	Dont need to run any of this
###############################
#some plots of hnet
library(ggplot2)

dat <- read.delim("~/HPC/Benchmarks/datasets/humanNet/humannet+thresh2.txt", sep="\t", header=TRUE, stringsAsFactors =FALSE)

x <- sort(dat$ms_uniprot_ac, decreasing=TRUE)
x = as.data.frame(cbind(x, labs= row.names(x)))
x = as.data.frame(table(x))
x[,1] = as.character(x[,1])
names(x) = c("uniprot_id", "Freq")
x = x[order(x$Freq, decreasing=TRUE),]
plot(x$Freq)


qplot(Freq, data=x[1:500,], geom="histogram", binwidth=100)  + theme(axis.text.x = element_text(angle=45))
 ggplot(x) + geom_histogram(aes(x=one, y=..count..))
 





















###############################
#Formatting TIP49 data


dat <- read.delim('~/Dropbox/kroganlab/dataset/apms_benchmark/datasets/TIP49/TIP49.txt', sep='\t', header=TRUE, stringsAsFactors =FALSE)

dat2 <- read.delim('~/Dropbox/kroganlab/dataset/apms_benchmark/datasets/TIP49/input_TIP49.txt', sep='\t', header=TRUE, stringsAsFactors =FALSE)


names(dat2) <- dat2[1,]
dat2 <- dat2[-c(1:2),]
dat2 <- dat2[,-c(2:4)]

#function that takes a graph/adjacency matrix and returns the long form of it
adjmat2long <- function(le_data){	
	tmp = cbind(Bait=le_data[,1], score = le_data[,2], Prey = names(le_data)[2])
	names(tmp) <- c("Bait", "score", "PREY")
	for(i in 3:(dim(le_data)[2])){
		tmp = rbind(tmp, cbind(Bait=le_data[,1], score = le_data[,i], Prey = names(le_data)[i]))
	}
	tmp = as.data.frame(tmp)
	tmp$score = as.numeric(tmp$score)
	return(tmp)
}

#get the long form of the TIP49 input data
dat2 <- adjmat2long(dat2)


















