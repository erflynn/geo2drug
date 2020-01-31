# get the lists of data where we have manual annotations

# // TODO: I can only seem to find this for human

require('tidyverse')
source("code/utils/labeling_utils.R")

manual <- read.delim("../data/manual/geo_manual_labels_jdw.tsv")
table(manual$gender)
sex_labeled_manual <- manual %>% 
  filter(gender %in% c("F", "M", "MIXED")) %>% 
  select(X, gender) %>% 
  rename(gsm="X", sex="gender")

gsm_filt_w_gse <- read.csv("data/01_sample_lists/gse_gsm_all_geo_dedup.csv")

# most overlap, some are missing from de-duplication -- I am ok w this
length(intersect(sex_labeled_manual$gsm, gsm_filt_w_gse$gsm)) 
length(setdiff(sex_labeled_manual$gsm, gsm_filt_w_gse$gsm))

manual_w_gse <- left_join(sex_labeled_manual, gsm_filt_w_gse)
manual_w_gse %>% select(gse) %>% unique() %>% dplyr::count() # 363 studies

# get platform data
con <- dbConnect(SQLite(), GEOMETADB_PATH)
gse.list <- formattedList(unique(manual_w_gse$gse))
gpls <- dbGetQuery(con, sprintf("SELECT gse, gpl FROM gse_gpl WHERE gse IN ('%s');", gse.list))
dbDisconnect(con)

gpls2 <-gpls %>% group_by(gse) %>% summarize(gpl=paste(gpl, collapse="|")) %>% ungroup()
manual_w_gpl <- left_join(manual_w_gse, gpls2, by="gse") %>% dplyr::rename(Gender=sex)

# now summarize to study
manual_gpl_summ <- summarizeToStudy(manual_w_gpl)

manual_gpl_summ %>% 
  filter(!study_type=="no_labels") %>% 
  filter(!is.na(gse)) %>%
  write_csv("data/01_sample_lists/human_manual_studies.csv")

# also grab the mixed sample! this is interesting...
manual_w_gpl %>% filter(!Gender %in% c("F", "M")) %>% dplyr::rename(sex=Gender) %>%
  filter(!is.na(gse)) %>% write_csv("data/01_sample_lists/human_mixed_sample_manual.csv")
