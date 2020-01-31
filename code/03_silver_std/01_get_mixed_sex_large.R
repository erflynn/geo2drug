# get a list of mixed sex large studies for silver-standard
require('tidyverse')
options(stringsAsFactors=FALSE)

human_gse_counts <- read.csv("data/01_sample_lists/human_gse_counts.csv")
table(human_gse_counts$study_type)

gse_to_keep <- read.delim("data/02_labeled_data/non_cell_line_gse.txt")
gse_no_cell <- filter(human_gse_counts, 
                      gse %in% gse_to_keep$gse_to_keep) %>% 
  filter(study_type == "both" & num_f >= 10 & num_m >= 10)

# --> 860 studies
write.csv(select(gse_no_cell, -study_type), file="data/01_sample_lists/gse_for_silver_std_human.csv", row.names=FALSE)

# check on the number of gsms in these data
gsm_filt_w_gse <- read.csv("data/01_sample_lists/gse_gsm_all_geo_dedup.csv")
gse_gsm <- inner_join(gse_no_cell, gsm_filt_w_gse[,c("gse", "gsm")])
length(unique(gse_gsm$gsm))

# mouse + rat
mouse_gse_counts <- read.csv("data/01_sample_lists/mouse_gse_counts.csv")
gse_mouse <- filter(mouse_gse_counts, gse %in% gse_to_keep$gse_to_keep) %>% filter(study_type == "both" & num_f >= 10 & num_m >= 10)
write.csv(select(gse_mouse, -study_type), file="data/01_sample_lists/gse_for_silver_std_mouse.csv", row.names=FALSE)

rat_gse_counts <- read.csv("data/01_sample_lists/rat_gse_counts.csv")
gse_rat <- filter(rat_gse_counts, gse %in% gse_to_keep$gse_to_keep) %>% filter(study_type == "both" & num_f >= 10 & num_m >= 10)
write.csv(select(gse_rat, -study_type), file="data/01_sample_lists/gse_for_silver_std_rat.csv", row.names=FALSE)
