
require('GEOmetadb')
require('tidyverse')
require('lubridate')
require('exprsex')

GEOMETADB_PATH <- "../GEOmetadb.sqlite" # downloaded 9/23/19

LIST_ORGANISMS <- list("human"="Homo sapiens", "mouse"="Mus musculus", "rat"="Rattus norvegicus")

#'Formats a list of items so it can be used for a SQL query
#'
#' @param my_list the list of items
#' @return formatted list, properly escaped and separated by quotes
formattedList <- function(my_list){
  paste(unique(my_list), collapse="\', \'")
}

reformatAleDat <- function(dat){
  dat2 <- dat %>% rename("gsm"="X", "gse"="ExperimentID", "gpl"="PlatformID") %>% 
    mutate(gsm=paste("GSM", gsm, sep=""), gse=paste("GSE", gse, sep=""), gpl=paste("GPL", gpl, sep=""))
  return(dat2)
}
summarizeToStudy <- function(dat){
  study_dat <- dat %>% group_by(gse) %>% 
    summarise(num_f=sum(Gender=="F"), num_m=sum(Gender=="M"), total=n(), gpl=paste(unique(gpl), collapse="|")) %>% 
    mutate(study_type = ifelse(num_f == 0, ifelse(num_m ==0, "no_labels", "m"), ifelse(num_m==0, "f", "both")))
  return(study_dat)
}