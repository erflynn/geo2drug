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