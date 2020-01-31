require('tidyverse')
require('GEOmetadb')
options(stringsAsFactors=FALSE)

human_gse_counts <- read.csv("data/01_sample_lists/human_gse_counts.csv")

gse_to_keep <- read.delim("data/02_labeled_data/non_cell_line_gse.txt")
f_only_no_cell <- filter(human_gse_counts, gse %in% gse_to_keep$gse_to_keep) %>% 
  filter(study_type == "f" & num_f >= 10 & num_f==total)
m_only_no_cell <- filter(human_gse_counts, gse %in% gse_to_keep$gse_to_keep) %>% 
  filter(study_type == "m" & num_m >= 10 & num_m==total)

# for now -- grab the same platforms as the train data
human_ss <- read.csv("data/01_sample_lists/gse_for_silver_std_human.csv")
gpl_counts <- human_ss %>% 
  separate_rows(gpl, sep="\\|") %>% 
  group_by(gpl) %>% 
  summarize(n=length(gse)) %>% 
  arrange(desc(n))

top_ten_gpl <- (gpl_counts %>% head(10) %>% select(gpl))$gpl
f_only <- sepReformatGPL(f_only_no_cell) %>% filter(gpl %in% top_ten_gpl)
m_only <- sepReformatGPL(m_only_no_cell) %>% filter(gpl %in% top_ten_gpl)

set.seed(2)

# grap up to 4 from each of the top ten platforms
f_only2 <- sampWrapper(f_only, 4)
m_only2 <- sampWrapper(m_only, 4)

rbind(f_only2, m_only2) %>% 
  select(gse, gpl, study_type) %>% 
  write_csv("data/01_sample_lists/human_single_sex_studies.csv")

