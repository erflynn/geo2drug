# get a list of mixed sex large studies for silver-standard
require('tidyverse')
options(stringsAsFactors=FALSE)


human_gse_counts <- read.csv("data/01_sample_lists/human_gse_counts.csv")
table(human_gse_counts$study_type)

gse_to_keep <- read.delim("data/02_labeled_data/non_cell_line_gse.txt")
gse_no_cell <- filter(human_gse_counts, 
                      gse %in% gse_to_keep$gse_to_keep) %>% 
  filter(study_type == "both" & num_f >= 10 & num_m >= 10)

seq_gpls_df <- read_csv("data/01_sample_lists/seq_gpls.csv")
gpls_seq_list <- c(seq_gpls_df$mouse, seq_gpls_df$human)



h_w_plat <- gse_no_cell %>% mutate(plat_type=ifelse(gpl %in% gpls_seq_list, "seq", "oligo"))
h_w_plat %>% group_by(plat_type) %>% count() # 1017 oligo, 71 seq
# --> 1088 studies
write.csv(h_w_plat %>% filter(plat_type=="oligo") %>% select(-study_type, -plat_type),  
          file="data/01_sample_lists/gse_for_silver_std_human_oligo.csv", row.names=FALSE)

write.csv(h_w_plat %>% filter(plat_type=="seq") %>% select(-study_type, -plat_type),  
          file="data/01_sample_lists/gse_for_silver_std_human_seq.csv", row.names=FALSE)

# check on the number of gsms in these data
gsm_filt_w_gse <- read.csv("data/01_sample_lists/gse_gsm_all_dedup.csv")
gse_gsm <- inner_join(gse_no_cell, gsm_filt_w_gse[,c("gse", "gsm")])
length(unique(gse_gsm$gsm))

# mouse 
mouse_gse_counts <- read.csv("data/01_sample_lists/mouse_gse_counts.csv")
gse_mouse <- filter(mouse_gse_counts, gse %in% gse_to_keep$gse_to_keep) %>%
  filter(study_type == "both" & num_f >= 10 & num_m >= 10) # 222

m_w_plat <- gse_mouse %>% mutate(plat_type=ifelse(gpl %in% gpls_seq_list, "seq", "oligo"))
m_w_plat %>% group_by(plat_type) %>% count() # 142 oligo, 80 seq
write.csv(m_w_plat %>% filter(plat_type=="oligo") %>% select(-study_type, -plat_type), 
          file="data/01_sample_lists/gse_for_silver_std_mouse_oligo.csv", row.names=FALSE)
write.csv(m_w_plat %>% filter(plat_type=="seq") %>% select( -study_type, -plat_type), 
          file="data/01_sample_lists/gse_for_silver_std_mouse_seq.csv", row.names=FALSE)
