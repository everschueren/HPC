library(sqldf)

vpu_mg132_ifn_evidence <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-IFN-MG132-Ub/052013-tlj-63-64-reip-evidence.txt", header=T)
vpu_evidence <- read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-Ub/input/051313-tlj-52-55-evidence.txt", header=T)

vpu_gfp_l = sqldf("select Modified_Sequence, Proteins, avg(log(Intensity_L)) as 'GFP' from vpu_evidence where Raw_file = 'VE20130513-04' or  Raw_file = 'VE20130513-05' group by Proteins")
vpu_gfp_h = sqldf("select Modified_Sequence, Proteins, avg(log(Intensity_H)) as 'GFP' from vpu_evidence where Raw_file = 'VE20130513-07' or  Raw_file = 'VE20130513-08' group by Proteins")
vpu_gfp = rbind(vpu_gfp_l, vpu_gfp_h)

vpu_vpu_h = sqldf("select Modified_Sequence, Proteins, avg(log(Intensity_H)) as 'VPU' from vpu_evidence where Raw_file = 'VE20130513-04' or  Raw_file = 'VE20130513-05' group by Proteins")
vpu_vpu_l = sqldf("select Modified_Sequence, Proteins, avg(log(Intensity_L)) as 'VPU' from vpu_evidence where Raw_file = 'VE20130513-07' or  Raw_file = 'VE20130513-08' group by Proteins")
vpu_vpu = rbind(vpu_vpu_l, vpu_vpu_h)

vpu_mg132_ifn_gfp_l = sqldf("select Modified_Sequence, Proteins, avg(log(Intensity_L)) as 'GFP' from vpu_mg132_ifn_evidence where Raw_file = 'VE20130520-37' or  Raw_file = 'VE20130520-38' group by Proteins")
vpu_mg132_ifn_gfp_h = sqldf("select Modified_Sequence, Proteins, avg(log(Intensity_H)) as 'GFP' from vpu_mg132_ifn_evidence where Raw_file = 'VE20130520-34' or  Raw_file = 'VE20130520-35' group by Proteins")
vpu_mg132_ifn_gfp = rbind(vpu_mg132_ifn_gfp_l, vpu_mg132_ifn_gfp_h)

vpu_mg132_ifn_vpu_h = sqldf("select Modified_Sequence, Proteins, avg(log(Intensity_H)) as 'VPU' from vpu_mg132_ifn_evidence where Raw_file = 'VE20130520-37' or  Raw_file = 'VE20130520-38' group by Proteins")
vpu_mg132_ifn_vpu_l = sqldf("select Modified_Sequence, Proteins, avg(log(Intensity_L)) as 'VPU' from vpu_mg132_ifn_evidence where Raw_file = 'VE20130520-34' or  Raw_file = 'VE20130520-35' group by Proteins")
vpu_mg132_ifn_vpu = rbind(vpu_mg132_ifn_vpu_l, vpu_mg132_ifn_vpu_h)

gfp_all_combined = sqldf("select G1.GFP as 'G1', G2.GFP as 'G2' from vpu_gfp G1 join vpu_mg132_ifn_gfp G2 on G1.Proteins = G2.Proteins ")
plot(gfp_all_combined$G1, gfp_all_combined$G2)

vpu_all_combined = sqldf("select V1.VPU as 'V1', V2.VPU as 'V2' from vpu_vpu V1 join vpu_mg132_ifn_vpu V2 on V1.Proteins = V2.Proteins")
plot(vpu_all_combined$V1, vpu_all_combined$V2)
