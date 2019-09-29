
require('GEOmetadb')
require('tidyverse')
require('lubridate')
require('exprsex')

GEOMETADB_PATH <- "../GEOmetadb.sqlite" # downloaded 9/23/19

LIST_ORGANISMS <- list("human"="Homo sapiens", "mouse"="Mus musculus", "rat"="Rattus norvegicus")
MIN_ROWS <- 10000 # minimum number of rows for the data

#'Formats a list of items so it can be used for a SQL query
#'
#' @param my_list the list of items
#' @return formatted list, properly escaped and separated by quotes
formattedList <- function(my_list){
  paste(unique(my_list), collapse="\', \'")
}
