library(RMySQL)
library(stringr)

## make db connection and get all data
## find a solution for the ssh tunnel that can't be run in R 
#ssh -L 3307:localhost:3306 everschueren@bluemoon.ucsf.edu

connect = function(){
  con = dbConnect(MySQL(),host="127.0.0.1",dbname="HPCKrogan",user="everschueren",pass="hpckr0gan",port=3307)  
  con
}

disconnect = function(con){
  dbDisconnect(con)
}

LoadProspectorData = function(attributeList="*",selectCondition=""){
  con = connect()
  query = str_join("select ", attributeList ," from MS_Result_prospector R join MS_Experiment E on E.ms_experiment_id = R.ms_experiment_id ", selectCondition)
  result = dbGetQuery(con,query)
  disconnect(con)
  result
}

