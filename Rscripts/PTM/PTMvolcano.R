vif_diff <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_reversed_labels/042313-tlj-59-62-differential.txt")
plot(vif_diff$VIF_GFP_logFC, -log10(vif_diff$VIF_GFP_adjPVal) , pch=20, ylab="-log10(FDR)", xlab = "log2 FC")
abline(h=-log10(0.2),col="darkgrey", lty=2, lwd=2)
abline(v=c(-1.5,1.5), col="darkgrey", lty=2, lwd=2)
