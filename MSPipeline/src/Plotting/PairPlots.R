#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(optparse))
suppressMessages(library(compiler))

PairPlot.panel.cor = function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method="pearson",use="pairwise.complete.obs"),2)
  txt <- format(c(r, 0.12), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  text(0.5, 0.5, txt, cex=1.5)
}

PairPlot.panel.points = function(x,y, ...){
#   d = cbind(x,y)
#   d = d[d[,1]!=0 & d[,2]!=0,]
  points(x,y,pch=20)
  f = lm(y~x)
  abline(f, col="red", lwd=1)
}

PairPlot.main = function(matrix_file, output_file, data_idx=3){
  all = read.delim(matrix_file, stringsAsFactors=F)
  data = all[,data_idx:ncol(all)]
  cnames = colnames(data)
  cnames = gsub("_","\n",cnames) ## replace all underscores by dashes
  cnames = gsub("(\\\n)(L|H)","_\\2",cnames) ## re-replace all \n followed by L or H with _L _H
  colnames(data) = cnames
  dmin = quantile(data, .0001 ,na.rm=T)[[1]]
  dmax = quantile(data, .9999 ,na.rm=T)[[1]]
  dlim = max(abs(dmin), abs(dmax))
  pdf(file=output_file, width=7, height=7)
  ## for variable axis 
  #pairs(data, upper.panel=PairPlot.panel.cor, lower.panel=PairPlot.panel.points, xlim=c(-dlim,dlim), ylim=c(-dlim,dlim), las=1, cex.axis=.6, cex.labels=1)
  ## fixed axis
  pairs(data, upper.panel=PairPlot.panel.cor, lower.panel=PairPlot.panel.points, xlim=c(-5,5), ylim=c(-5,5), las=1, cex.axis=.7, cex.labels=1)
  dev.off()
}

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-d", "--data_file"),
              help="data file containing value matrix"),
  make_option(c("-o", "--output_file"),
              help="output file for pair plot PDF")
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
PairPlot.main(matrix_file=parsedArgs$data_file, output_file=parsedArgs$output_file)  

#PairPlot.main(matrix_file="~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-all-Ub/processed_repeats/Mock_v_WT_all_evidence_FLT_MAT.txt", output_file="~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Infection-PTM/Mock-v-WT-all-Ub/processed_repeats/Mock_v_WT_all_evidence_COR.pdf")

