library(limma)

silacRI.plot = function(ratios, L_intensities, H_intensities){
  VAL_idx = 3
  cols = c("blue", "red")
  linecols = c("darkblue", "darkred")
  for(i in VAL_idx:ncol(ratios)){
    R = ratios[, c(1:(VAL_idx-1),i)]
    R[, VAL_idx] = log2(R[,VAL_idx])
    I = L_intensities[, c(1:(VAL_idx-1),i)]
    I[,VAL_idx] = (log2(I[,VAL_idx]) + log2(H_intensities[, i]))/2 
  
    RI = merge(R, I, by=c("uniprot_id","Sequence"))
    RI.omitNA = na.omit(RI)
    #RI.omitNA = RI
    
    if(i == VAL_idx){
      fun=plot
    }else{
      fun=points
    }
    fun(RI.omitNA[,VAL_idx+1], RI.omitNA[,VAL_idx], pch=20, col=cols[i-VAL_idx+1], xlab="Avg(Log2(L)+Log2(H))", ylab="Log2(L/H)")
    
    I.fit = lowess(RI.omitNA[,VAL_idx+1], RI.omitNA[,VAL_idx], f=0.4)
    lines(I.fit$x, I.fit$y, lwd=3, col=linecols[i-VAL_idx+1])
  }
  legend("topright", legend=colnames(ratios)[VAL_idx:ncol(ratios)], fill=cols, cex=.6)
}

silacRI.normalize = function(df, method="scale"){
  VAL_idx = 3
  df[, VAL_idx:ncol(df)] = as.data.frame(normalizeBetweenArrays(as.matrix(df[, VAL_idx:ncol(df)], method=method)))
  df
}

# return z-scores
silacRI.IDZscore = function(ratios, H_intensities, L_intensities){
  VAL_idx = 3
  z_scores = ratios[,1:(VAL_idx-1)]
  for(i in VAL_idx:ncol(ratios)){
    col_idxs = c(1:(VAL_idx-1),i)
    R = ratios[, col_idxs]
    R = na.omit(R)
    R[, VAL_idx] = log2(R[,VAL_idx])
    
    IL = L_intensities[, col_idxs]
    IL = na.omit(IL)
    IH = H_intensities[, col_idxs]
    I[,VAL_idx] = (log2(I[,VAL_idx]) + log2(H_intensities[, i]))/2 
    
    RI = merge(R, I, by=c("uniprot_id","Sequence"), all=T)
    RI_omitNA = na.omit(RI)
    colnames(RI_omitNA)[3:4] = c("Ratio", "Intensity") 
    RI_omitNA = sqldf("select * from RI_omitNA order by Intensity ASC")
    
    bin_size = 200
    bins=ceiling(nrow(RI_omitNA)/bin_size) ## 200 is a good bin size
    break_points = c(1,(seq(1,bins) * bin_size))
    break_points[length(break_points)] = nrow(RI_omitNA)
    break_values = RI_omitNA$Intensity[break_points]
    
    RI_omitNA = cbind(RI_omitNA, cut(RI_omitNA$Intensity, break_values, include.lowest=T,right=F))
    colnames(RI_omitNA)[5] = "bin" 
    RI_omitNA_bins = sqldf("select bin, avg(Ratio) as 'I_mean', stdev(Ratio) as 'I_sd' from RI_omitNA group by bin")
    RI_2 = sqldf("select RI.uniprot_id, RI.Sequence, RI.Ratio / RIB.I_sd as 'IDZ' from RI_omitNA RI left join RI_omitNA_bins RIB on RI.bin= RIB.bin")
    z_scores = merge(z_scores, RI_2, by=c("uniprot_id", "Sequence"), all=T) 
  }
  
  colnames(z_scores) = colnames(ratios)
  z_scores
}

ratios = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_ratio_intensity/042313-tlj-59-62-matrix-ratios-ubi_unique_filtered.txt")
plot(log2(ratios$Vif_H_GFP_L), log2(ratios$GFP_H_Vif_L), pch=20, col="gray", xlim=c(-4,4), ylim=c(-4,4), xlab="Log2(Vif_H / GFP_L)", ylab="Log2(GFP_H / Vif_L)")
abline(h=0,lwd=2, col="darkgrey")
abline(v=0,lwd=2, col="darkgrey")
L_intensities = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_ratio_intensity/042313-tlj-59-62-matrix-L_intensities-ubi_unique_filtered.txt")
H_intensities = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_ratio_intensity/042313-tlj-59-62-matrix-H_intensities-ubi_unique_filtered.txt")


pdf("~/Projects/HPCKrogan/Data/HIV-proteomics/Plots/Ub-scans/Jurkat_VIF_RI_plot.pdf", width=7, height=5)
silacRI.plot(ratios, L_intensities, H_intensities)
dev.off()

ratios = silacRI.normalize(ratios)
L_intensities = silacRI.normalize(L_intensities)
H_intensities = silacRI.normalize(H_intensities)
z_scores = silacRI.IDZscore(ratios, L_intensities, H_intensities)
write.table(z_scores, file="../../Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_ratio_intensity/042313-tlj-59-62-matrix-z_scores.txt", eol="\n", sep="\t", quote=F, row.names=F, col.names=T)


pdf("~/Projects/HPCKrogan/Data/HIV-proteomics/Plots/Ub-scans/Jurkat_VIF_RI_normalized_plot.pdf", width=7, height=5)
silacRI.plot(ratios, L_intensities, H_intensities)
dev.off()

ratios = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpr-Ub/processed_ratio_intensity/042313-tlj-57-60-matrix-ratios-ubi_unique_filtered.txt")
L_intensities = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpr-Ub/processed_ratio_intensity/042313-tlj-57-60-matrix-L_intensities-ubi_unique_filtered.txt")
H_intensities = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpr-Ub/processed_ratio_intensity/042313-tlj-57-60-matrix-H_intensities-ubi_unique_filtered.txt")

pdf("~/Projects/HPCKrogan/Data/HIV-proteomics/Plots/Ub-scans/Jurkat_VPR_RI_plot.pdf", width=7, height=5)
silacRI.plot(ratios, L_intensities, H_intensities)
dev.off()

ratios = silacRI.normalize(ratios)
L_intensities = silacRI.normalize(L_intensities)
H_intensities = silacRI.normalize(H_intensities)

pdf("~/Projects/HPCKrogan/Data/HIV-proteomics/Plots/Ub-scans/Jurkat_VPR_RI_normalized_plot.pdf", width=7, height=5)
silacRI.plot(ratios, L_intensities, H_intensities)
dev.off()

ratios = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-Ub/processed_ratio_intensity/051313-tlj-52-55-matrix-ratios-ubi_unique_filtered.txt")
L_intensities = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-Ub/processed_ratio_intensity/051313-tlj-52-55-matrix-L_intensities-ubi_unique_filtered.txt")
H_intensities = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-Ub/processed_ratio_intensity/051313-tlj-52-55-matrix-H_intensities-ubi_unique_filtered.txt")

pdf("~/Projects/HPCKrogan/Data/HIV-proteomics/Plots/Ub-scans/Jurkat_VPU_RI_plot.pdf", width=7, height=5)
silacRI.plot(ratios, L_intensities, H_intensities)
dev.off()

ratios = silacRI.normalize(ratios)
L_intensities = silacRI.normalize(L_intensities)
H_intensities = silacRI.normalize(H_intensities)

pdf("~/Projects/HPCKrogan/Data/HIV-proteomics/Plots/Ub-scans/Jurkat_VPU_RI_plot_normalized.pdf", width=7, height=5)
silacRI.plot(ratios, L_intensities, H_intensities)
dev.off()
