# annot_geo_metadata.R
# E Flynn
# Last Updated: 4/16/2019
#
# Grab all the GEOMetadb data for mouse, rat, and human at the GSM, GSE level.
# Deduplicate by mapping each sample to the earliest study it is present in.

require('GEOmetadb')
require('lubridate')
require('tidyverse')

con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite')

# download all MOUSE, RAT, HUMAN METADATA at the gsm level
gsms <- dbGetQuery(con, 'SELECT * FROM gsm ;')

gsm_filt <- filter(gsms,
  organism_ch1 %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus"),
  molecule_ch1== "total RNA",
  channel_count==1)

gsm_filt <- gsm_filt[,c("gsm", "series_id", "gpl", "title", "description", "source_name_ch1", "characteristics_ch1")] %>% 
  rename("gse"="series_id")
gsm_filt2 <- separate_rows(gsm_filt, gse, sep=",")

# grab the study data for these gsms 
formatted_list_gses <- paste(unique(gsm_filt2$gse), collapse="\', \'")
gses <- dbGetQuery(con, sprintf("SELECT * FROM gse WHERE gse IN ('%s');", formatted_list_gses))
gse_filt <- gses[,c("gse", "title", "summary", "pubmed_id", "submission_date", "overall_design")]
dbDisconnect(con)

# map a sample to study with the earliest date if it is duplicated
sample_study_date <- inner_join(gsm_filt2[,c("gsm", "gse")], select(gse_filt, c("gse", "submission_date")))
sample_study_date$submission_date <- as.Date(sample_study_date$submission_date)
mapped_to_earliest <- sample_study_date %>% group_by(gsm) %>% summarize(gse= gse[which.min(submission_date)])  # SUPER SLOW
write.csv(mapped_to_earliest, file="data/db_data/gse_gsm_all_geo_dedup.csv", row.names=FALSE)

# write out lists of all the gsm + gse metadata 
gsm_filt_w_gse <- inner_join(select(gsm_filt, -gse), mapped_to_earliest, by="gsm")
gse_data <- gse_filt[gse_filt$gse %in% mapped_to_earliest$gse,]
write.csv(gse_data, file="data/db_data/gse_all_geo_info.csv", row.names=FALSE)
write.csv(gsm_filt_w_gse, file="data/db_data/gsm_all_geo_info.csv", row.names=FALSE)



