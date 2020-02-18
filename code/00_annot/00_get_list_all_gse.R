# get_list_all_gse.R
# E Flynn
# Last Updated: 9/23/2019
#
# Grab all the GEOMetadb data for mouse and human at the GSM, GSE level.
# Deduplicate by mapping each sample to the earliest study it is present in.

source("code/utils/labeling_utils.R")


# download all data at the gsm level and filter for containing the species of interest

con <- dbConnect(SQLite(), GEOMETADB_PATH)

gsms <- dbGetQuery(con, 'SELECT * FROM gsm ;') %>% 
  select("gsm", "series_id", "gpl", "organism_ch1", "title", "description", "source_name_ch1", "characteristics_ch1", "molecule_ch1", "channel_count") %>%
  rename("gse"="series_id") %>%
  filter(str_detect(organism_ch1, "Homo sapiens|Mus musculus"),
         molecule_ch1 == "total RNA",
         channel_count == 1) %>% 
  select(-channel_count, -molecule_ch1) %>%
  separate_rows(organism_ch1, sep=";") %>%
  as_tibble() 

dbDisconnect(con)


# clean up the names
# - for many of the Mus musculus names it specifies the subspecies --> we are going to just rename these to Mus musculus
# - remove tabs and spaces
gsms2 <- gsms %>%
  mutate(organism_ch1 = ifelse(str_detect(organism_ch1, "Mus musculus"), 
                            "Mus musculus", organism_ch1)) %>% 
  mutate(organism_ch1 = ifelse(str_detect(organism_ch1, "Homo sapiens"), 
                             "Homo sapiens", organism_ch1))  %>%
  filter(organism_ch1 %in% c("Homo sapiens", "Mus musculus")) 


# ---  divide by species --- #
MIN.ROWS <- 10000
# NOTE: by filtering for at least 10k rows, we remove all the HTS data
LIST.OLIGO.TECH <- c("in situ oligonucleotide", "spotted oligonucleotide", "spotted DNA/cDNA", "oligonucleotide beads")
# this removes MPSS, RT-PCR


con <- dbConnect(SQLite(), GEOMETADB_PATH)

human <- filter(gsms2, organism_ch1=="Homo sapiens") %>% 
  select(gse, gsm, everything())
human_gpls <- dbGetQuery(con, sprintf("SELECT gpl, organism, technology, data_row_count FROM gpl WHERE gpl IN ('%s');",
                                      formattedList(unique(human$gpl)))) %>% 
  unique() # 2001 
human_gpls2 <- human_gpls %>% filter(data_row_count >= MIN.ROWS)  %>%
  mutate(organism = ifelse(str_detect(organism, "Homo sapiens"), 
                           "Homo sapiens", organism)) # 1005
human_gpls3 <- human_gpls2 %>% filter(organism=="Homo sapiens") # 969
human_gpls4 <- human_gpls3 %>% filter(technology %in% LIST.OLIGO.TECH) # 963
write_csv(human_gpls4, "data/01_sample_lists/human_gpl.csv")
human_gsms <- human %>% filter(gpl %in% human_gpls4$gpl) # 754131
human_gsms %>% write_csv("data/01_sample_lists/human_gsms.csv")

human_gses <- human_gsms %>% separate_rows(gse, sep=",") %>% select(gse) %>% unique() # 19837

mouse <- filter(gsms2, organism_ch1=="Mus musculus") %>% 
  select(gse, gsm, everything())
mouse_gpls <- dbGetQuery(con, sprintf("SELECT gpl, organism, technology, data_row_count FROM gpl WHERE gpl IN ('%s');",
                                      formattedList(unique(mouse$gpl)))) %>% 
  unique() # 1035 
dbDisconnect(con)

mouse_gpls2 <- mouse_gpls %>% filter(data_row_count >= MIN.ROWS)  %>%
  mutate(organism = ifelse(str_detect(organism, "Mus musculus"), 
                           "Mus musculus", organism)) # 481
mouse_gpls3 <- mouse_gpls2 %>% filter(organism=="Mus musculus") # 451
mouse_gpls4 <- mouse_gpls3 %>% filter(technology %in% LIST.OLIGO.TECH) # 440

write_csv(mouse_gpls4, "data/01_sample_lists/mouse_gpl.csv")
mouse_gsms <- mouse %>% filter(gpl %in% mouse_gpls4$gpl) # 216599

mouse_gsms %>% write_csv("data/01_sample_lists/mouse_gsms.csv")
mouse_gses <- mouse_gsms %>% separate_rows(gse, sep=",") %>% select(gse) %>% unique() # 14797

# there is a tiny bit of overlap - 22 samples
human_mouse <- intersect(human_gsms$gsm, mouse_gsms$gsm) 
write_csv(data.frame(human_mouse),"data/01_sample_lists/multiple_organism_data.csv" )

oligo_gse_list <- c(human_gses$gse, mouse_gses$gse) 


# look at the metadata for HTS studies -- this is from ARCHS4
miceadds::load.Rdata("data/00_db_data/human_gsm_meta.rda", "human_meta_seq")
miceadds::load.Rdata("data/00_db_data/mouse_gsm_meta.rda", "mouse_meta_seq")

human_seq_gsms <- names(human_meta_seq) # 238,522
mouse_seq_gsms <- names(mouse_meta_seq) # 284,907


# if there are multiple GSEs for a GSM, they are in a list --> this works
human_seq_gses <- unique(unlist(lapply(human_meta_seq, function(x) x$Sample_series_id))) # 9124
mouse_seq_gses <- unique(unlist(lapply(mouse_meta_seq, function(x) x$Sample_series_id))) # 9642

human_seq_gpls <- unique(unlist(lapply(human_meta_seq, function(x) x$Sample_platform_id))) # 8
mouse_seq_gpls <- unique(unlist(lapply(mouse_meta_seq, function(x) x$Sample_platform_id))) # 8
cbind("human"=human_seq_gpls, "mouse"=mouse_seq_gpls) %>% data.frame() %>% write_csv("data/01_sample_lists/seq_gpls.csv")



seq_gse_list <- c(human_seq_gses, mouse_seq_gses)
gse_list <- c(oligo_gse_list, seq_gse_list) # 53400

# ----grab the study data for these gsms ---- #
# we want this for BOTH seq and microarray data
con <- dbConnect(SQLite(), GEOMETADB_PATH)

gses <- dbGetQuery(con, sprintf("SELECT * FROM gse WHERE gse IN ('%s');", 
                                formattedList(gse_list)))
gse_filt <- gses %>% 
  select("gse", "title", "pubmed_id", "submission_date", "overall_design",  "summary") 
# 49686

# add a GPL column
gse_gpl <- dbGetQuery(con, sprintf("SELECT gse, gpl FROM gse_gpl WHERE gse IN ('%s');", 
                                   formattedList(gse_filt$gse)))
dbDisconnect(con)
gse_gpl2 <- gse_gpl %>% group_by(gse) %>% dplyr::summarize(gpl=paste(gpl, collapse=",")) %>% ungroup()
gse_filt2 <- left_join(gse_filt, gse_gpl2)

gse_filt3 <- gse_filt2 %>% mutate(study_type=case_when(
  gpl %in% human_seq_gpls ~ "human_seq",
  gpl %in% mouse_seq_gpls ~ "mouse_seq",
  gpl %in% human_gpls4 ~ "human_oligo",
  gpl %in% mouse_gpls4 ~ "mouse_oligo"
))
gse_filt4 <- gse_filt3 %>% separate(study_type, into=c("organism", "plat_type"), sep="_")

gse_filt5 <- gse_filt4 %>% select(gse, gpl, organism, plat_type, everything())
write_csv(gse_filt5, "data/01_sample_lists/gse_metadata_all.csv")

# -----  get duplicated and non-duplicated GSMs ----- #

micro <- rbind((mouse_gsms %>% select(gse, gsm)), (human_gsms %>% select(gse, gsm)))
micro2 <- micro %>% separate_rows(gse, sep=",") 
con <- dbConnect(SQLite(), GEOMETADB_PATH)
seq_gse_gsm <- dbGetQuery(con, sprintf("SELECT gse, gsm FROM gse_gsm WHERE gsm IN ('%s');", 
                                   formattedList(c(human_seq_gsms, mouse_seq_gsms))))

dbDisconnect(con)

# map a sample to study with the earliest date if it is duplicated
sample_study_date <- inner_join(select(rbind(micro2, seq_gse_gsm), c("gsm", "gse")),
                                select(gse_filt4, c("gse", "submission_date"))) %>%
  mutate(submission_date=ymd(submission_date))

mapped_to_earliest <- sample_study_date %>%
  group_by(gsm) %>%
  dplyr::summarize(gse= gse[which.min(submission_date)]) # SLOW

write_csv(mapped_to_earliest, "data/01_sample_lists/gse_gsm_all_dedup.csv")

# now do a little filtering
gse_list2 <- unique(mapped_to_earliest$gse) # 41941 / 53400

gse_filt5 <- gse_filt4 %>% filter(gse %in% gse_list2)
write_csv(gse_filt5, "data/01_sample_lists/gse_metadata_all_filt.csv")


