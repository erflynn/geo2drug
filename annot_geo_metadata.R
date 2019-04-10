# annot_geo_metadata.R
# 
# Grab all the GEOMetadb data 
require('GEOmetadb')
require('lubridate')
require('tidyverse')


con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite')

# download all MOUSE, RAT, HUMAN METADATA --> somewhere
gsms <- dbGetQuery(con, 'SELECT * FROM gsm ;')

gsm_filt <- filter(gsms,
  organism_ch1 %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus"),
  molecule_ch1== "total RNA",
  channel_count==1)

# assign each to the earlier study that contains it?

gsm_filt <- gsm_filt[,c("gsm", "series_id", "gpl", "title", "description", "source_name_ch1", "characteristics_ch1")] %>% rename("gse"="series_id")
gsm_filt2 <- separate_rows(gsm_filt, gse, sep=",")

# grab the study data
formatted_list_gses <- paste(unique(gsm_filt2$gse), collapse="\', \'")
gses <- dbGetQuery(con, sprintf("SELECT * FROM gse WHERE gse IN ('%s');", formatted_list_gses))
gse_filt <- gses[,c("gse", "title", "summary", "pubmed_id", "submission_date", "overall_design")]
dbDisconnect(con)

# map to the earliest date!
sample_study_date <- inner_join(gsm_filt2[,c("gsm", "gse")], select(gse_filt, c("gse", "submission_date")))
sample_study_date$submission_date <- as.Date(sample_study_date$submission_date)
mapped_to_earliest <- sample_study_date %>% group_by(gsm) %>% summarize(gse= gse[which.min(submission_date)])  # SUPER SLOW
# TODO: only map duplicated!                            
write.csv(mapped_to_earliest, file="data/db_data/gse_gsm_all_geo_dedup.csv", row.names=FALSE)

gsm_filt_w_gse <- inner_join(select(gsm_filt, -gse), mapped_to_earliest, by="gsm")
gse_data <- gse_filt[gse_filt$gse %in% mapped_to_earliest$gse,]
write.csv(gse_data, file="data/db_data/gse_all_geo_info.csv", row.names=FALSE)
write.csv(gsm_filt_w_gse, file="data/db_data/gsm_all_geo_info.csv", row.names=FALSE)

# separate into study vs sample
gse_level_info <- gsm_filt_w_gse %>% group_by(gse) %>% summarize(study_description=length(unique(description))==1,
                                                                      study_source_name_ch1=length(unique(source_name_ch1))==1, 
                                                                      study_characteristics_ch1=length(unique(characteristics_ch1))==1 )

gsm_info_w_keep <- inner_join(gsm_filt_w_gse, gse_level_info,by="gse")

sample_strs <- lapply(1:nrow(gsm_info_w_keep),
       function(i) {
        paste(gsm_info_w_keep[i, "title"], paste( gsm_info_w_keep[i,c("description", "source_name_ch1", "characteristics_ch1")
                   ][!gsm_info_w_keep[i,c("study_description", "study_source_name_ch1", "study_characteristics_ch1")]], 
      collapse=" | "), sep=" | ")})

## do this for study level
gse_level_info[,c("study_description", "study_source_name_ch1", "study_characteristics_ch1")]

