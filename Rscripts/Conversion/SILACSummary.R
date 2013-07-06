library(reshape2)

input = read.delim("~/Projects/HPCKrogan/Data/HIV-proteomics/Jurkat-Expression-PTM/Vif-Ub/input/042313-tlj-59-62-evidence.txt")
input = input[,c("Modified.sequence","Raw.file","Ratio.H.L")]
input[is.na(input)]=1

input_reshaped = dcast(input, Modified.sequence ~ Raw.file, value.var="Ratio.H.L", fun.aggregate=mean)
input_reshaped[is.na(input_reshaped)]=1

boxplot(log2(input_reshaped[,2:5]), las =2)

plot(density(log2(input_reshaped[,2]))) 
for(i in 3:5){
  lines(density(log2(input_reshaped[,i]))) 
}
peptides_num = length(unique(input$Sequence))
