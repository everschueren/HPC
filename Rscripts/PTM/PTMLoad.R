library(RMySQL)

## make db connection and get all data
## find a solution for the ssh tunnel that can't be run in R 
#system("ssh -L 3307:localhost:3306 everschueren@bluemoon.ucsf.edu")
con = dbConnect(MySQL(),host="127.0.0.1",dbname="HPCKrogan",user="everschueren",pass="hpckr0gan",port=3307)

MS_Experiment = dbGetQuery(con, "select * from MS_Experiment")
MS_Result_simple = dbGetQuery(con, "select * from MS_Result_maxquant")
MS_Result_combined = dbGetQuery(con, "select * from MS_Result_maxquant R join MS_Experiment E on E.ms_experiment_id = R.ms_experiment_id")
  
#MS_Result_simple_pairs = dbGetQuery(con, "select * from MS_Result_simple R1 join MS_Result_simple R2 on R1.peptide=R2.peptide join MS_Experiment E1 on R1.ms_experiment_id=E1.id join MS_Experiment E2 on R2.ms_experiment_id=E2.id where R1.peptide like '%(gl)%' and R1.sample != R2.sample and R1.silac_heavy = R2.silac_heavy and R1.silac_light = R2.silac_light and E1.cell_overexpressed_proteins = 20261 and E2.cell_overexpressed_proteins = 20261")

dbDisconnect(con)