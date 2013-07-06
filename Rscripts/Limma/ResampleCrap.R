



















# # compute average score and counts
# value_name = colnames(input)[VAL_idx]
# input_grouped = sqldf(str_join("select ",id_name,", avg(",value_name,") as 'average_score', count(*) as 'count' from input group by ",id_name," order by average_score desc"))
# 
# P_abundance_1 = resample.Peptides.C(scored_data=input_grouped, id_idx=1, score_idx=2, count_idx=3, sample_size=SAMPLE_SIZE)
# 
# VAL_idx = 4
# value_name = colnames(input)[VAL_idx]
# input_grouped$average_score = input_grouped$average_score * -1 
# input_grouped = sqldf(str_join("select ",id_name,", avg(",value_name,") as 'average_score', count(*) as 'count' from input group by ",id_name," order by average_score desc"))
# P_abundance_2 = resample.Peptides.C(scored_data=input_grouped, id_idx=1, score_idx=2, count_idx=3, sample_size=SAMPLE_SIZE)
# 
# input_grouped_acc = cbind(input_grouped, P_abundance_1, P_abundance_2)
# plot(input_grouped_acc$P_abundance_1, input_grouped_acc$P_abundance_2)
# 
# 
# 
# 
# 
# 
# 
# input_values_grouped_selected = sqldf(str_join("select * from input_grouped_acc where average_score > 0 and P_abundance_1 <= ",PVALUE," or P_abundance_2 <= ",PVALUE," order by P_abundance_1, P_abundance_2 asc"))
# input_values_grouped_selected_annotated = annotate_with_uniprot(data=input_values_grouped_selected,key="uniprot_id","HUMAN")
# 
# 
# 
# ## visualize_scores as volcano plot
# plot(input_values_grouped$average_score, -log2(input_values_grouped$P_abundance), pch=20)
# abline(h=-log2(PVALUE), col = "lightgray", lty = 3, lwd=3)
# 
# # input = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vpu-Ub/processed/051313-tlj-52-55-differential.txt")
# #input = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/processed_reversed_labels/042313-tlj-59-62-differential.txt")
# # input_values = input[,c(ID_idx,LFC_idx)]
# # value_name = colnames(input)[LFC_idx]
# # # compute average score and counts 
# # input_values_grouped = sqldf(str_join("select uniprot_ac, avg(",value_name,") as 'average_score', count(*) as 'count' from input_values group by uniprot_ac order by average_score desc"))
# # P_abundance = resample.Peptides.C(scored_data=input_values_grouped, id_idx=1, score_idx=2, count_idx=3, sample_size=SAMPLE_SIZE)
# # input_values_grouped = cbind(input_values_grouped, P_abundance)
# # input_values_grouped_selected = sqldf(str_join("select * from input_values_grouped where average_score > 0 and P_abundance <= ",PVALUE," order by P_abundance asc, average_score desc"))
# # 
# # ## visualize_scores as volcano plot
# # plot(input_values_grouped$average_score, -log2(input_values_grouped$P_abundance), pch=20)
# # abline(h=-log2(PVALUE), col = "lightgray", lty = 3, lwd=3)
