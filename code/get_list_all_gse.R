# get_list_all_gse.R
# E Flynn
# Last Updated: 9/23/2019
#
# Grab all the GEOMetadb data for mouse, rat, and human at the GSM, GSE level.
# Deduplicate by mapping each sample to the earliest study it is present in.

source("code/utils/labeling_utils.R")


# download all data at the gsm level and filter for containing the species of interest

con <- dbConnect(SQLite(), GEOMETADB_PATH)

gsms <- dbGetQuery(con, 'SELECT * FROM gsm ;') %>% 
  select("gsm", "series_id", "gpl", "organism_ch1", "title", "description", "source_name_ch1", "characteristics_ch1", "molecule_ch1", "channel_count") %>%
  rename("gse"="series_id") %>%
  filter(str_detect(organism_ch1, "Homo sapiens|Mus musculus|Rattus norvegicus"),
         molecule_ch1 == "total RNA",
         channel_count == 1) %>% 
  select(-channel_count, -molecule_ch1) %>%
  separate_rows(organism_ch1, sep=";") %>%
  as_tibble() 

# is the RNA-seq GPL in here?
# if it is - keep the metadata for these

# data_row_count > 10k | data_row_count==0
# do we want to filter out other technology?
#

dbDisconnect(con)


# clean up the names
# - for many of the Mus musculus names it specifies the subspecies --> we are going to just rename these to Mus musculus
# - remove tabs and spaces
gsms2 <- gsms %>%
  mutate(organism_ch1 = ifelse(str_detect(organism_ch1, "Mus musculus"), 
                            "Mus musculus", organism_ch1)) %>% 
  mutate(organism_ch1 = ifelse(str_detect(organism_ch1, "Rattus norvegicus"), 
                                 "Rattus norvegicus", organism_ch1)) %>%
  mutate(organism_ch1 = ifelse(str_detect(organism_ch1, "Homo sapiens"), 
                             "Homo sapiens", organism_ch1))  %>%
  filter(organism_ch1 %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus")) %>%
  separate_rows(gse, sep=",") # separate out by study


# ---  divide by species --- #



human <- filter(gsms2, organism_ch1=="Homo sapiens") %>% 
  select(gse, gsm, everything())
rat <- filter(gsms2, organism_ch1=="Rattus norvegicus") %>% 
  select(gse, gsm, everything())
mouse <- filter(gsms2, organism_ch1=="Mus musculus") %>% 
  select(gse, gsm, everything())

# --- get GPL platform mapping --- #

#' Get the organisms for a list of GPLs
#'
#' @param df - the dataframe containing a gpl column
#' @param org_name - the organism name
#' @param con - connection to GEOMetadb
#' 
#' @return map_gpls - the org mapping to GPLs
getOrgMatchingGPL <- function(df, org_name, con){
  org_gpl <- df %>% select(gpl) %>% unique()
  map_gpls <- dbGetQuery(con, sprintf("SELECT gpl, organism, data_row_count FROM gpl WHERE gpl IN ('%s');",
                                      formattedList(org_gpl$gpl))) %>%
    unique() %>% 
    filter(data_row_count >= MIN_ROWS) %>%
    mutate(organism = ifelse(str_detect(organism, org_name),
                                 org_name, organism)) %>%
    as_tibble() 
  matching_dat <- inner_join(map_gpls, df, by="gpl") %>%
    filter(organism_ch1==organism) %>%
    unique()
  return(matching_dat)
}
con <- dbConnect(SQLite(), GEOMETADB_PATH)

hs <- getOrgMatchingGPL(human, "Homo sapiens", con)
rn <- getOrgMatchingGPL(rat, "Rattus norvegicus", con)
mm <- getOrgMatchingGPL(mouse, "Mus musculus", con)

dbDisconnect(con)

human_rat <- intersect(hs$gsm, rn$gsm) # none
human_mouse <- intersect(hs$gsm, mm$gsm)
rat_mouse <- intersect(mm$gsm, rn$gsm) # none

# mult_org_dat
#   gsm gse gpl organism mult_resolved?
write_csv(data.frame(human_mouse),"data/sample_lists/multiple_organism_data.csv" )

gse_list <- c(rn$gse, hs$gse, mm$gse) 

# ----grab the study data for these gsms ---- #

con <- dbConnect(SQLite(), GEOMETADB_PATH)

gses <- dbGetQuery(con, sprintf("SELECT * FROM gse WHERE gse IN ('%s');", 
                                formattedList(gse_list)))
gse_filt <- gses %>% 
  select("gse", "title", "summary", "pubmed_id", "submission_date", "overall_design")

dbDisconnect(con)

# map a sample to study with the earliest date if it is duplicated
sample_study_date <- inner_join(select(gsms2, c("gsm", "gse")), 
                                select(gse_filt, c("gse", "submission_date"))) %>% 
  mutate(submission_date=ymd(submission_date))

mapped_to_earliest <- sample_study_date %>% 
  group_by(gsm) %>% 
  summarize(gse= gse[which.min(submission_date)]) # SLOW

write_csv(mapped_to_earliest, "data/sample_lists/gse_gsm_all_geo_dedup.csv")


gse_list2 <- unique(mapped_to_earliest$gse)


# write out the gse metadata
gse_data <- gses %>% 
  select("gse", "title", "pubmed_id", "submission_date", "overall_design",  "summary") %>% 
  filter(gse %in% gse_list2)

# add in an organism column?
combined_gsm_gse <- do.call(rbind, list(mm, rn, hs))
gse_data_org <- left_join(gse_data, select(combined_gsm_gse, gse, organism_ch1), by="gse") %>%
  select(gse, organism_ch1, everything()) %>%
  rename(organism=organism_ch1) %>%
  unique()  %>%
  as_tibble()


write_csv(gse_data_org, "data/sample_lists/gse_all_geo_info.csv")

# filter and write out the mouse, rat, and human data
mm2 <- inner_join(mm, mapped_to_earliest, by=c("gsm", "gse"))  %>% 
  select(-organism, -organism_ch1, -data_row_count) %>% 
  select( gse, gsm, everything())
rn2 <- inner_join(rn, mapped_to_earliest, by=c("gsm", "gse")) %>% 
  select(-organism, -organism_ch1, -data_row_count) %>% 
  select( gse, gsm, everything())
hs2 <- inner_join(hs, mapped_to_earliest, by=c("gsm", "gse")) %>% 
  select(-organism, -organism_ch1, -data_row_count) %>% 
  select( gse, gsm, everything())


mm2 %>% write_csv("data/sample_lists/mouse_gse_gsm.csv")
rn2 %>% write_csv("data/sample_lists/rat_gse_gsm.csv")
hs2 %>% write_csv("data/sample_lists/human_gse_gsm.csv")

mm2 %>% select(gse) %>% unique() %>% write_csv("data/sample_lists/gse_mouse.csv")
rn2 %>% select(gse) %>% unique() %>% write_csv("data/sample_lists/gse_rat.csv")
hs2 %>% select(gse) %>% unique() %>% write_csv("data/sample_lists/gse_human.csv")


# write out a list of all the GPLs
con <- dbConnect(SQLite(), GEOMETADB_PATH)

mm_gpls <- dbGetQuery(con, sprintf("SELECT * FROM gpl WHERE gpl IN ('%s');",
                                    formattedList(mm2$gpl)))
rn_gpls <- dbGetQuery(con, sprintf("SELECT * FROM gpl WHERE gpl IN ('%s');",
                                   formattedList(rn2$gpl)))
hs_gpls <- dbGetQuery(con, sprintf("SELECT * FROM gpl WHERE gpl IN ('%s');",
                                   formattedList(hs2$gpl)))
dbDisconnect(con)

hs_gpls %>% 
  select(gpl, title, manufacturer, technology, description, data_row_count) %>% 
  write_csv("data/sample_lists/human_gpl.csv")

rn_gpls %>% 
  select(gpl, title, manufacturer, technology, description, data_row_count) %>% 
  write_csv("data/sample_lists/rat_gpl.csv")

mm_gpls %>% 
  select(gpl, title, manufacturer, technology, description, data_row_count) %>% 
  write_csv("data/sample_lists/mouse_gpl.csv")


# which are the most common GPLs for each?

# what is the consensus gene list for these?
# what are the top gpls for each?

mm_count_gpls <- mm2 %>% select(gse, gpl) %>% unique() %>% group_by(gpl) %>% mutate(n=n()) %>% select(-gse) %>% unique() %>% arrange(desc(n))
rn_count_gpls <- rn2 %>% select(gse, gpl) %>% unique() %>% group_by(gpl) %>% mutate(n=n()) %>% select(-gse) %>% unique() %>% arrange(desc(n))
hs_count_gpls <- hs2 %>% select(gse, gpl) %>% unique() %>% group_by(gpl) %>% mutate(n=n()) %>% select(-gse) %>% unique() %>% arrange(desc(n))

# get the genes for each of these
mm_count_gpls$gpl[1:8]



## ---- SCRATCH CODE ---- ##

rnaseq_gpls <- dbGetQuery(con, "SELECT * FROM gpl WHERE organism=='Homo sapiens' AND technology=='high-throughput sequencing';")
#rnaseq_gpls <- dbGetQuery(con, "SELECT * FROM gpl WHERE organism=='Homo sapiens' AND technology=='high-throughput sequencing';")

query_alt <- dbGetQuery(con, 'SELECT * FROM gsm ;') %>% 
  select("gsm", "series_id", "gpl", "organism_ch1", "title", "description", "source_name_ch1", "characteristics_ch1", "molecule_ch1", "channel_count") %>%
  rename("gse"="series_id") %>%
  filter(str_detect(organism_ch1, "Homo sapiens|Mus musculus|Rattus norvegicus")) %>%
  as_tibble()

# what happens when we change the filtering
gpls <- query_alt %>% select(gpl) %>% unique()
gpl_info <- dbGetQuery(con, sprintf("SELECT gpl, organism, title, manufacturer, technology, description, data_row_count FROM gpl WHERE gpl IN ('%s');",
                                    formattedList(gpls$gpl))) %>%
  filter(data_row_count >= 10000,
         str_detect(organism, "Homo sapiens|Mus musculus|Rattus norvegicus")
         ) %>%
  as_tibble()

gpl_info %>% group_by(technology) %>% count()

# ---- start with platforms --- #
all_gpls <- dbGetQuery(con, sprintf("SELECT * from gpl;"))
all_gpls %>% 
  filter(data_row_count >= 10000 | data_row_count==0) %>% 
  group_by(technology) %>% count()


gsms <- dbGetQuery(con, 'SELECT * FROM gsm WHERE gpl=="GPL15433";') %>% 
  select("gsm", "series_id", "gpl", "organism_ch1", "title", "description", "source_name_ch1", "characteristics_ch1", "molecule_ch1", "channel_count") %>%
  rename("gse"="series_id") 
