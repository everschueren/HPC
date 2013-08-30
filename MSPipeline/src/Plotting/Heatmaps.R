suppressMessages(library(gplots))
suppressMessages(library(SDMTools))

heatmap.EV = function(data, cexRow=1.0, palette.breaks=NULL, color.panel=NULL){
  
  if(is.null(color.panel)){
    color.panel = colorpanel(10, "gold", "white", "steelblue")  
  }
  if(is.null(palette.breaks)){
    palette.breaks = seq(from=0,to=1,by=.1)  
  }
  
  margins = c(1,1)
  lmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
  lwid = c(0.5,4.5,2)
  lhei = c(0.5,4.5,2)
  
  #heatmap.2(data, dendrogram="row",key=F,trace="none",col=color.panel,labCol=rep("",ncol(data)),lmat=lmat, lhei=lhei, lwid=lwid, margins=margins, cexRow=cexRow, breaks=palette.breaks)
  
  heatmap.2(data, dendrogram="row",key=F,trace="none",col=color.panel,lmat=lmat, lhei=lhei, lwid=lwid, margins=margins, cexRow=cexRow, breaks=palette.breaks, cexCol=cexRow)
  
  ## print legend
  leg_x_offset = 0.95
  leg_y_offset = 0.01
  legend_coords = cbind(x = c(0,0.015,0.015,0) + leg_x_offset, y =c(0.1,0.1,0,0)  + leg_y_offset)
  legend.gradient(legend_coords,cols=color.panel, limits=c(round(min(palette.breaks),2),round(max(palette.breaks),2)), title="",cex=cexRow)
}
