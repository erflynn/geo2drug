# Code for identifying which mat files do not convert

require('tidyverse')
args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]

# // TODO - do not need GSE-GPLS for step 00, need this for step 01
silver_std <- read.csv(sprintf("data/01_sample_lists/gse_for_silver_std_%s.csv", organism))
silver_std_s <- silver_std %>% filter(!str_detect(gpl, "\\|")) 
silver_std_m <- silver_std %>% filter(str_detect(gpl, "\\|"))  %>% separate_rows(gpl, sep="\\|") %>%
  mutate(gse=sprintf("%s-%s", gse, gpl)) 

silver_std_gse <- c(silver_std_s$gse, silver_std_m$gse)

downloaded_f <- list.files(sprintf("data/03_silver_std/%s/00_mat_files/", organism))
gse.gpls <- unique(sapply(downloaded_f, function(x) strsplit(x, "_")[[1]][[1]])) # w GPL suffix
missing.gse.gpls <- setdiff(silver_std_gse, gse.gpls)
data.frame(missing.gse.gpls) %>% write_csv(sprintf("data/03_silver_std/%s/mat_missing.csv", organism))

gses <- unique(sapply(downloaded_f, function(x) strsplit(x, "_|-")[[1]][[1]])) # w/o GPL suffix
print(sprintf("For %s, successfully downloaded %s gses, %s gse-gpls out of %s in the silver std (gses-gpls=%s)", 
              organism, length(gses), length(gse.gpls), length(silver_std$gse), length(silver_std_gse)))

# // TODO - expand to do this across steps
