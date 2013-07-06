idx_vector = c(9:16,41:48)
idx_vector = c(73:88)
test = bigMatrix[,idx_vector]
# colnames(test) = rep(c("gfp_h","vifmg132_h","vifmg132_l","gfp_l","gfp_h","vprmg132_h","vprmg132_l","gfp_l"),each=2)
colnames(test) = rep(c("vpu1_h","gfp_h","gfp_l","vpu1_l","vpu2_h","gfp_h","gfp_l","vpu2_l"),each=2)
biolrep=rep(1:8,each=2)
design <- model.matrix(~ 0+factor(rep(c(1,2,2,1,3,2,2,3),each=2)))
colnames(design) = c("vpu1","gfp","vpu2")
contrasts <- makeContrasts("(vpu1+vpu2)-(gfp)",levels=design)
result = doLimma(test, design, biolrep, contrasts)
result = sqldf("select * from result where logFC > 0")
result = cbind(1:nrow(result),result)
colnames(result)[1] = "rank"
result_annotated = sqldf(str_join("select * from result R left join peptide_protein_descriptions P on P.peptide = R.ID group by P.peptide, P.key order by R.P_Value asc"))

