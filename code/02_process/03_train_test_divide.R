# divide the data into training and testing

require('tidyverse')

args <- commandArgs(trailingOnly = TRUE)
organism <- args[1]

# get a list of all the files
list2 <- list.files(sprintf("data/03_silver_std/%s/02_keep_labels/", organism))

# read in the files
all_dat <- lapply(1:length(list2), function(i){
  f1 <- read.csv(sprintf("data/03_silver_std/%s/02_keep_labels/%s", organism, list2[[i]]))
  if(nrow(f1) < 10){ # deal with the empty files
    return(NA)
  }
  f1$gse <- strsplit(list2[[i]], "_")[[1]][[1]]
  return(f1)
})


all_dat2 <- do.call(rbind, all_dat[!is.na(all_dat)]) # 170 studies, 18460 samples

mult_gpl <- all_dat2 %>% 
  filter(str_detect(gse, "-")) %>% 
  select(gse) %>% 
  unique() %>% 
  separate(gse, "-", into=c("gse", "gpl"))
single_gpl <- all_dat2 %>%
  filter(!str_detect(gse, "-")) %>% select(gse) %>% unique()

gse_info <- read_csv(sprintf("data/01_sample_lists/gse_for_silver_std_%s.csv", organism))

gse_info_s <- gse_info %>% filter(gse %in% single_gpl$gse)
# add these -- they are in the path and we can't have them not separated out!
gse_info_m <- gse_info %>% separate_rows(gpl, sep="\\|") %>% semi_join(mult_gpl)
gse_info_m2 <- gse_info_m %>% mutate(gse=sprintf("%s-%s", gse, gpl))

gse_info2 <- rbind(gse_info_s, gse_info_m2)

set.seed(2)
gse_info2$rand.int <- sample(1:nrow(gse_info2), nrow(gse_info2), replace=FALSE)
training <- gse_info2 %>% group_by(gpl) %>% top_n(n=2, wt=rand.int) 
testing <- gse_info2 %>% filter(!gse %in% training$gse)
print(length(unique(testing$gse)))
print(length(unique(training$gse)))

all_dat2 %>% 
  filter(gse %in% training$gse) %>%
  mutate(gse=unlist(gse)) %>%
  select(gse, gsm, text_sex) %>% 
  write_csv(sprintf("data/01_sample_lists/%s_training.csv", organism))

all_dat2 %>% 
  filter(gse %in% testing$gse) %>% 
  mutate(gse=unlist(gse)) %>%
  select(gse, gsm, text_sex) %>% 
  write_csv(sprintf("data/01_sample_lists/%s_testing.csv", organism))
