require('tidyverse')
require('GEOmetadb')
options(stringsAsFactors=FALSE)

# // TODO:
#  - do this for mouse too

args <- commandArgs(trailingOnly = TRUE)
organism <- args[1]

gse_counts <- read.csv(sprintf("data/01_sample_lists/%s_gse_counts.csv", organism))
seq_gpls <- read_csv("data/01_sample_lists/seq_gpls.csv") %>% select(organism)
# bin and count
gse_counts %>% mutate(size_grp=ifelse(total > 1000, ">1000", ifelse(total > 100, ">100", ifelse(total > 10, ">10", "0-10")))) %>%
  group_by(study_type, size_grp) %>% mutate("organism"=organism) %>% count() %>% 
  write_csv(sprintf("data/01_sample_lists/%s_gse_counts_size.csv", organism))

gse_to_keep <- read.delim("data/02_labeled_data/non_cell_line_gse.txt")
f_only_no_cell <- filter(gse_counts, gse %in% gse_to_keep$gse_to_keep) %>% 
  filter(study_type == "f" & num_f >= 10 & num_f==total)
m_only_no_cell <- filter(gse_counts, gse %in% gse_to_keep$gse_to_keep) %>% 
  filter(study_type == "m" & num_m >= 10 & num_m==total)

f_w_plat <- f_only_no_cell %>% mutate(plat_type=ifelse(gpl %in% unlist(seq_gpls[,1]), "seq", "oligo"))
f_w_plat %>% group_by(plat_type) %>% count() # 490 oligo, 56 seq for human; 596 oligo, 75 seq for mouse

m_w_plat <- m_only_no_cell %>% mutate(plat_type=ifelse(gpl %in% unlist(seq_gpls[,1]), "seq", "oligo"))
m_w_plat %>% group_by(plat_type) %>% count() # 457 oligo, 82 seq for human; 1092 oligo, 150 seq for mouse

ss_plat <- rbind(f_w_plat, m_w_plat) 
ss_plat %>% filter(plat_type=="oligo") %>% write_csv(sprintf("data/01_sample_lists/%s_oligo_single_sex.csv", organism))
ss_plat %>% filter(plat_type=="seq") %>% write_csv(sprintf("data/01_sample_lists/%s_seq_single_sex.csv", organism))



##### -------------- STOP ----------------- #####



# for now -- grab the same platforms as the train data
ss <- read.csv(sprintf("data/01_sample_lists/gse_for_silver_std_%s.csv", organism))
gpl_counts <- ss %>% 
  separate_rows(gpl, sep="\\|") %>% 
  group_by(gpl) %>% 
  summarize(n=length(gse)) %>% 
  arrange(desc(n))

top_ten_gpl <- (gpl_counts %>% head(10) %>% select(gpl))$gpl
gpl_counts_ref %>% filter(gpl %in% top_ten_gpl)
setdiff(gpl_counts_ref[1:15,]$gpl, top_ten_gpl)



f_only <- sepReformatGPL(f_only_no_cell) %>% filter(gpl %in% top_ten_gpl)
m_only <- sepReformatGPL(m_only_no_cell) %>% filter(gpl %in% top_ten_gpl)

set.seed(2)

# grap up to 4 from each of the top ten platforms
f_only2 <- sampWrapper(f_only, 4)
m_only2 <- sampWrapper(m_only, 4)

rbind(f_only2, m_only2) %>% 
  select(gse, gpl, study_type) %>% 
  write_csv("data/01_sample_lists/%s_single_sex_studies.csv")

rbind(f_only, m_only) %>% 
  select(gse, gpl, study_type) %>% 
  write_csv(sprintf("data/01_sample_lists/%s_single_sex_studies_all.csv", organism))

# now write these out at a sample level for processing 
list.gses <- c(f_only2$gse, m_only2$gse)
ale <- read_csv(sprintf("data/01_sample_lists/%s_ale_sex_lab.csv", organism))
ale2 <- ale %>% 
  rename(text_sex="Gender") %>%
  mutate(text_sex=ifelse(text_sex=="M", "male", "female"))

sp_ale <- ale2 %>% filter(gse %in% list.gses)
mp_ale <- ale2 %>% mutate(gse=paste(gse, gpl, sep="-")) %>% filter(gse %in% list.gses)

samp_lab <-rbind(sp_ale, mp_ale) %>% select(gse, gsm, text_sex, gpl)
samp_lab %>% write_csv(sprintf("data/01_sample_lists/%s_single_sex_samples.csv", organism))
