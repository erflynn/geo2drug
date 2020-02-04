# divide the data into training and testing

require('tidyverse')

args <- commandArgs(trailingOnly = TRUE)
organism <- args[1]

# look at what files exist
downloaded_f <- list.files(sprintf("data/03_silver_std/%s/02_keep_labels/", organism))
gse.gpls <- unique(sapply(downloaded_f, function(x) strsplit(x, "_")[[1]][[1]])) # w GPL suffix


# GPL breakdown
silver_std <- read.csv(sprintf("data/01_sample_lists/gse_for_silver_std_%s.csv", organism), stringsAsFactors = FALSE)
silver_std_s <- silver_std %>% filter(!str_detect(gpl, "\\|")) 
silver_std_m <- silver_std %>% filter(str_detect(gpl, "\\|"))  %>% separate_rows(gpl, sep="\\|") %>%
  mutate(gse=sprintf("%s-%s", gse, gpl)) 

silver_std_gse <- c(silver_std_s$gse, silver_std_m$gse)
missing.gse.gpls <- setdiff(silver_std_gse, gse.gpls)
gses <- unique(sapply(downloaded_f, function(x) strsplit(x, "_|-")[[1]][[1]])) # w/o GPL suffix

print(sprintf("For %s, successfully downloaded %s gses, %s gse-gpls out of %s in the silver std (gses-gpls=%s)", 
              organism, length(gses), length(gse.gpls), length(silver_std$gse), length(silver_std_gse)))


silver_std2 <- rbind(silver_std_s, silver_std_m)
silver_std_rem <- silver_std2 %>% filter(gse %in% gse.gpls) 
gpl_counts <- silver_std_rem %>% group_by(gpl) %>% count() %>% arrange(desc(n)) 


# --- CONSTRUCT THE OUTPUT --- #
list2 <- list.files(sprintf("data/03_silver_std/%s/02_keep_labels/", organism))

# read in the files
all_dat <- lapply(1:length(list2), function(i){
  f1 <- read.csv(sprintf("data/03_silver_std/%s/02_keep_labels/%s", organism, list2[[i]]), 
                 stringsAsFactors = FALSE)
  if(nrow(f1) < 10){ # deal with the empty files
    return(NA)
  }
  f1$gse <- strsplit(list2[[i]], "_")[[1]][[1]]
  return(f1)
})

all_dat2 <- do.call(rbind, all_dat[!is.na(all_dat)]) # 593 studies, 70481 samples

all_dat3 <- all_dat2 %>% left_join(silver_std2 %>% select(gse, gpl))

writeOutput <- function(my.gses, lab){
  all_dat3 %>% 
    filter(gse %in% my.gses) %>%
    mutate(gse=unlist(gse)) %>%
    select(gse, gsm, text_sex, gpl) %>% 
    write_csv(sprintf("data/01_sample_lists/%s_%s.csv", organism, lab))
}

set.seed(2)

#### ---- code for human ---- ####
if (organism=="human"){
  mult_gpls <- gpl_counts %>% filter(n >=4)
  # there are *23* platforms with >=4 GPLs
  
  top10_gpl <- mult_gpls %>% head(10)
  full_gpl <- setdiff(mult_gpls$gpl, top10_gpl$gpl)
  
  # common: top 10 platforms: 4 train + 4 test
  rand_sel <- silver_std_rem %>% 
    filter(gpl %in% top10_gpl$gpl) %>% 
    group_by(gpl) %>%
    sample_n(8)
  train_common <- (rand_sel %>% sample_n(n=4))$gse
  test_common <- setdiff(rand_sel$gse, train_common)
  writeOutput(train_common, "training")
  writeOutput(test_common, "testing")
  
  # add in the rest of the GPLs with >=4, 2-3 train, 2-4 test
  n2_gpls <- mult_gpls %>% filter(n < 8 & gpl %in% full_gpl)
  n4_gpls <- mult_gpls %>% filter(n >=8 & gpl %in% full_gpl)
  n4_d <- silver_std_rem %>%
    filter(gpl %in% n4_gpls$gpl) %>%
    group_by(gpl) %>%
    sample_n(7) %>% 
    ungroup()
  dat_full<- rbind(n4_d, silver_std_rem %>%
                     filter(gpl %in% n2_gpls$gpl))
  test_full <- (dat_full %>% group_by(gpl) %>% sample_frac(0.5))$gse
  train_full <- setdiff(dat_full$gse, test_full)
  writeOutput(train_full, "training_full")
  writeOutput(test_full, "testing_full")
  
  # the rest of the studies with >= 4
  extended_test_dat <- silver_std_rem %>% 
    filter(gpl %in% mult_gpls$gpl & 
             !(gse %in% c(train_full, test_full, train_common, test_common)))
  extended_test <- extended_test_dat$gse
  writeOutput(extended_test, "testing_extended")
  
  # platforms with less than 4 studies
  rare_gpls <- silver_std_rem %>%
    filter(gpl %in% (gpl_counts %>% filter(n < 4))$gpl)
  extended_few_samples <- rare_gpls$gse
  writeOutput(extended_few_samples, "testing_rare")
}

#### ---- code for rat ---- ####
if (organism == "rat"){
  test_frac <- silver_std_rem %>% 
    group_by(gpl) %>%
    sample_frac(0.5) 
  train_gse <- setdiff(silver_std_rem$gse, test_frac$gse)
  test_gse <- test_frac$gse
  writeOutput(train_gse, "training")
  writeOutput(test_gse, "testing")
}

#### ---- code for mouse ---- ####
if (organism=="mouse"){
  # 24 mouse platforms, 9 have more than 1
  # split the 9 that have more than 1
  n2_gpls <- (filter(gpl_counts, n >1))$gpl
  train_df <- silver_std_rem %>% 
    filter(gpl %in% n2_gpls) %>%
    group_by(gpl) %>%
    sample_frac(0.5)
  writeOutput(train_df$gse, "training")
  
  test_df <- silver_std_rem %>% 
    filter(gpl %in% n2_gpls) %>% 
    filter(!gse %in% train_df$gse)
  writeOutput(test_df$gse, "testing")
  
  # write all the ones that have 1 to test_rare
  ext_test <- silver_std_rem %>% 
    filter(gpl %in% (filter(gpl_counts, n==1))$gpl)
  writeOutput(ext_test$gse, "testing_rare")
  
}



