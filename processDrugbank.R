# processDrugbank.R
# E Flynn
# 3/27/2019
# 
# Code for processing drugbank data --> data frame (condensed format) + vocab (long format)

require('rjson')
require('tidyverse')
source('drug_count_utils.R')

# load drugbank data and convert it into a data frame with one row per drug
drugbank <- fromJSON(file="data/db_data/drugbank_info.json") 
list.drugs <- lapply(drugbank, function(x) unique(unlist(tolower(c(x$synonyms, x$name)))))
names(list.drugs) <- lapply(drugbank, function(y) y$dbID)

drugbank_df <- do.call(rbind,
                       lapply(drugbank, function(x) data.frame(
                         lapply(x[c("synonyms","unii","name","cas","dbID","chebi" )], 
                                function(y) paste(y, collapse=" | ")), stringsAsFactors=FALSE)))
# escape quotes
drugbank_df$name <- sapply(drugbank_df$name, function(x) gsub( '\\"', "\\\'\'", x))
drugbank_df$synonyms <- sapply(drugbank_df$synonyms, function(x) gsub( '\\"', "\\\'\'", x))
write.table(drugbank_df, file="data/db_data/drugbank_parsed.txt", sep="\t", row.names=FALSE)


# --- CONSTRUCT DRUGBANK VOCABULARY --- #
drug_name_syn <- createNameSynVocab(drugbank_df, "dbID")
write.table(drug_name_syn, file="data/db_data/drugbank_vocab.txt", sep="\t", row.names=FALSE)
