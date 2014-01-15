## CONFIG

LOOKBACK=5 ## number of IPs to look back for carryover

## from text
## summary = read.delim("~/Box Documents/Projects/Herpes/summary/20140113/KSHV_errors_carryover.txt", stringsAsFactors=F)

## from SQL
library(RMySQL)

## execute this to access remote server
## ssh -L 3307:localhost:3306 higgs.ucsf.edu

con <- dbConnect(MySQL(), user="everschueren", password="", dbname="KroganMS", host="127.0.0.1", port=3307)
q = "select M.id, M.`bait_name`, D.`ms_protein_name`, D.`ms_num_unique_peptide` from MT_data D join MT_keys_merged M on D.id=M.id where M.`ip_cell_type` like 
'293%' and M.`data_set` like '%HHV8_glaunsinger%' and D.`ms_protein_name` like '%HHV8%'"
summary = dbGetQuery(con, q)
dbDisconnect(con)

## cleanup prey names
summary$ms_protein_name = gsub("\t|\n","",summary$ms_protein_name)
## cleanup bait names
summary$bait_name = gsub("_MG132|_mutant","",summary$bait_name) 

unique_ids = unique(summary$id)
annotations = c()
prev_factors = c()
bait_errors = c()

for(i in 1:length(unique_ids)){
  ## subselect
  id = unique_ids[i]
  sample = summary[summary$id == id,]
  
  ## run through sample entries to scan whether bait is detected or whether entry is carryover 
  bait_flag=F
  for(e in 1:nrow(sample)){
    entry = sample[e,]
    
    ## bait detected test
    bait_detected = grepl(entry$bait_name, entry$ms_protein_name)
    if(bait_detected){
      annotation = "BAIT"
      prev_factor=NA
      bait_flag=T
    }else{
      annotation = "CARRYOVER"
      actual_level = entry$ms_num_unique_peptide
      
      ## look back for previous levels
      j = i-LOOKBACK
      if(j < 1) j=1
      
      max_prev_level = actual_level
        
      while(j < i){
        prev_id = unique_ids[j]
        prev_sample = summary[summary$id == prev_id,]
        for(f in 1:nrow(prev_sample)){
          prev_entry = prev_sample[f,]
          ## bait detected test
          prev_bait_detected = grepl(prev_entry$bait_name, entry$ms_protein_name)
          if(prev_bait_detected){
            prev_level = prev_entry$ms_num_unique_peptide
            if(prev_level > max_prev_level) max_prev_level=prev_level
          }
        }
        prev_factor = max_prev_level/actual_level
        j=j+1
      }
    }
    annotations = c(annotations, annotation)
    prev_factors = c(prev_factors, prev_factor)
  }
  if(bait_flag){
    bait_error = rep("OK" ,nrow(sample))
  }else{
    bait_error = rep("BAIT NOT FOUND" ,nrow(sample))
  }
  bait_errors = c(bait_errors, bait_error)  
}

summary = data.frame(summary, annotations=annotations, carryover_factors=prev_factors, ip_errors=bait_errors) 
write.table(summary, file="Box Documents/Projects/Herpes/summary/20140113/KSHV_errors_carryover_report.txt", eol="\n", sep="\t", quote=F, row.names=F, col.names=T)

