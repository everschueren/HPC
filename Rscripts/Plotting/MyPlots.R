library(stringr)
library(gplots)
library(SDMTools)

#######################################################
## FUNCTIONS ##########################################

## boxplot with room for labels
boxplot.EV = function(data, main="", cex=1.0, xlab="", ylab=ylab){
  op = par(no.readonly = TRUE)
  par(mar = c(9, 4, 4, 2) + 0.1)
  boxplot(data, main=main, varwidth=T, xaxt="n", range=1, pch=20, cex=.5, ylab=ylab, las=1)
  axis(1, labels = FALSE, at=1:ncol(data))
  labels = colnames(data)
  text(1:ncol(data), par("usr")[3] - 0.75, srt = 45, adj = 1, labels = labels, xpd = TRUE, cex=cex, las =1)
  mtext(1, text = xlab, line = 8)
  par = op
}

error.bar = function(x, y, upper, lower=upper, length=0.1){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)){
    stop("vectors must be same length")
  }
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length)
}

## barplot with room for labels
barplot.EV = function(data, main="", cex=1.0, xlab="", ylab="", errorvalues=0, ylim_max=0, col="lightgray", beside=T){
  if(ylim_max==0){
    ylim_max = max(data+errorvalues)  
  }
  op = par(no.readonly = TRUE)
  par(mar = c(9, 4, 4, 2) + 0.1)
  xses = barplot(as.matrix(data), main=main, axes=F, axisnames=F, las=1, ylab=ylab, ylim=c(0,ylim_max), col=col, beside=beside)
  if(is.matrix(errorvalues)==T){
    error.bar(xses, data, errorvalues)
  }
  labels = colnames(data)
  # print(labels)
  text(xses, par("usr")[3] - 0.75, srt = 45, adj = 1, labels = labels, xpd = TRUE, cex=cex)
  #mtext(1, text = xlab, line = 8)
  axis(2, las=2)
  par = op
}

heatmap.EV = function(data){
  yellowwhiteblue = colorpanel(20, "gold", "white", "steelblue")
  margins = c(1,1)
  lmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
  lwid = c(0.5,4,2.5)
  lhei = c(0.5,4,2.5) 
  heatmap.2(data, dendrogram="row",key=F,trace="none",col=yellowwhiteblue,labCol=rep("",ncol(data)),lmat=lmat, lhei=lhei, lwid=lwid, margins=margins)
  ## print legend
  leg_x_offset = 0.05
  leg_y_offset = 0.05
  legend_coords = cbind(x =c(0,0.05,0.05,0) + leg_x_offset, y =c(0.15,0.15,0,0)  + leg_y_offset)
  legend.gradient(legend_coords,cols=yellowwhiteblue, limits=c(round(min(data),2),round(max(data),2)), title="")
}


