# now get a filtered list by platforms / platform frequency
require('tidyverse')
source("code/utils/assessment_utils.R")
source("code/utils/labeling_utils.R")

human_gse_counts <- read_csv("data/01_sample_lists/human_gse_counts.csv") %>% 
  mutate(organism="human")
mouse_gse_counts <- read_csv("data/01_sample_lists/mouse_gse_counts.csv") %>% 
  mutate(organism="mouse")
gse_counts <- rbind(human_gse_counts, mouse_gse_counts)

seq_gpls_df <- read_csv("data/01_sample_lists/seq_gpls.csv")
gpls_seq_list <- c(seq_gpls_df$mouse, seq_gpls_df$human)

con <- dbConnect(SQLite(), GEOMETADB_PATH)
gpls <- dbGetQuery(con, sprintf("SELECT gpl,technology FROM gpl WHERE gpl IN ('%s');", 
                                formattedList(unique(gse_counts$gpl))))

dbConnect(con)

gse_counts2 <- gse_counts %>% 
  mutate(plat_type=ifelse(gpl %in% gpls_seq_list, "seq", "oligo"))

# why are there this many without mapping platforms?
counts_by_plat <- gse_counts2 %>% 
  group_by(organism, plat_type, gpl) %>% count() %>% 
  arrange(desc(n)) 

# all the SEQ platforms
# the top 10 oligo platforms for each organism

ttLabel <- function(ds){
  ns <- nrow(ds)
  num_train <- ceiling(ns/2)
  return(c(rep(0, num_train), rep(1, ns-num_train)))
}

# REP_SAMPLE FIRST
getStudyLists <- function(plat){
  plat_data <- gse_counts2 %>% 
    filter(gpl %in% plat)
  
  rep_sample <- plat_data %>% 
    filter(study_type!="no_labels") %>%
    sample_n(10) %>%
    mutate(grp="rep_sample") %>%
    select(gse, grp, gpl)
  rep_sample$tt <- ttLabel(rep_sample)
  
  rem_data <- plat_data %>% 
    filter(!gse %in% rep_sample$gse)
  
  # Mixed sex + single sex 2nd
  mixed_sex <- rem_data %>% 
    filter(num_f >= 10 & num_m >= 10) %>%
    sampWrapper(12) %>%
    mutate(grp="mixed_sex") %>%
    select(gse, grp, gpl)
  mixed_sex$tt <- ttLabel(mixed_sex)
  
  single_sex_f <- rem_data %>% 
    filter(num_f >= 10 & total==num_f) %>%
    sampWrapper(8) %>%
    mutate(grp="ss_f") %>%
    select(gse, grp, gpl)
  single_sex_f$tt <- ttLabel(single_sex_f)
  
  single_sex_m <- rem_data %>%
    filter(num_m >= 10 & total==num_m) %>%
    sampWrapper(8) %>%
    mutate(grp="ss_m") %>%
    select(gse, grp, gpl)
  single_sex_m$tt <- ttLabel(single_sex_m)
  
  my_df <- do.call(rbind, list(rep_sample, mixed_sex, 
                               single_sex_f, single_sex_m))

  return(my_df)
}

plat.human <- counts_by_plat %>% 
  filter(plat_type=="oligo" & organism=="human") %>%
  filter(n >= 8) %>% ungroup() %>% select(gpl) 

plat.mouse <-counts_by_plat %>% 
  filter(plat_type=="oligo" & organism=="mouse") %>% 
  filter(n >= 8) %>% ungroup() %>% select(gpl)


top.plat.human <- plat.human %>% head(10)
top.plat.mouse <- plat.mouse %>% head(10)

set.seed(4)
h.lists <- lapply(top.plat.human$gpl, function(x) getStudyLists(c(x)))
m.lists <- lapply(top.plat.mouse$gpl, function(x) getStudyLists(c(x)))

h.df <- do.call(rbind, h.lists) 
h.df %>% write_csv("data/01_sample_lists/human_oligo_grps.csv")
m.df <- do.call(rbind, m.lists) 
m.df %>% write_csv("data/01_sample_lists/mouse_oligo_grps.csv")


# create train/valid/test + rep_sample for seq but not divided by sample
h.seq <- counts_by_plat %>% 
  filter(plat_type=="seq" & organism=="human") %>% 
  ungroup()
m.seq <- counts_by_plat %>% 
  filter(plat_type=="seq" & organism=="mouse") %>% 
  ungroup()

h.seq.lists <- getStudyLists(h.seq$gpl)
m.seq.lists <- getStudyLists(m.seq$gpl)
write_csv(h.seq.lists, "data/01_sample_lists/human_seq_grps.csv")
write_csv(m.seq.lists, "data/01_sample_lists/mouse_seq_grps.csv")


### get a representative sample for other platforms

repSampleOnly <- function(plat){
  plat_data <- gse_counts2 %>% 
    filter(gpl==plat)
  
  rep_sample <- plat_data %>% 
    filter(study_type!="no_labels") %>%
    sampWrapper(4) %>%
    mutate(grp="rep_sample") %>%
    select(gse, grp, gpl)
  return(rep_sample)
}

set.seed(6)
h.other.plat.rep <- lapply(setdiff(plat.human$gpl,top.plat.human$gpl), function(x) repSampleOnly(c(x)))
h.other <- do.call(rbind, h.other.plat.rep) %>% 
  write_csv("data/01_sample_lists/human_other_plat_rep.csv")
m.other.plat.rep <- lapply(setdiff(plat.mouse$gpl,top.plat.mouse$gpl), function(x) repSampleOnly(c(x)))
m.other <- do.call(rbind, m.other.plat.rep) %>% 
  write_csv("data/01_sample_lists/mouse_other_plat_rep.csv")

