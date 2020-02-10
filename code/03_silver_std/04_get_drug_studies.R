# code for getting drug studies from CREEDS, BROAD + our annotations
#   result:  organism, GSE, drug, annot_source
options(stringsAsFactors = FALSE)

# ---------- look at the NER results -------------- #
ner_annot <- read_csv("data/02_labeled_data/02_drug_labeling/ner_annot.csv")
ner_annot$X1 <- NULL
drug_mentions <-ner_annot %>% filter(chemical != "set()") %>% mutate(chemical=str_replace_all(chemical, "\\{|\\}",  "")) %>% separate_rows(chemical, sep=",") 
drug_mentions2 <- 
  drug_mentions %>% mutate(chemical=tolower(str_trim(chemical))) %>%
  mutate(chemical=str_replace_all(chemical, "'", "")) %>%
  filter((!str_detect(chemical, "gsm|gse"))) %>%
  filter(nchar(chemical) >= 3)

chemicals <-unique(drug_mentions2$chemical)

overlap <- intersect((drug_info_df %>% filter(numchar>3))$name,chemicals)
# 1219 overlap... I'm surprised?
novel <- setdiff(chemicals, (drug_info_df %>% filter(numchar>3))$name )
# many of these *ARE* drugs
#  but we grab them b/c we tokenize differently
#  HOWEVER; there are issues with misspellings/abbreviations + special characters (e.g. alpha)
# worth looking into these more

# -------------------------------------------------- #

# our results
drugbank_mapped <- read.delim2("data/02_labeled_data/drugbank_mapped_gse.txt", sep="\t")
# combine creeds, BROAD 
drug_pert_auto <- read.csv("data/00_db_data/single_drug_perturbations-p1.0.csv")
creeds_auto <- drug_pert_auto %>% filter(organism !="rat") %>% as_tibble()

drug_pert_manual <- read.csv("data/00_db_data/single_drug_perturbations-v1.0.csv")
creeds_manual <- drug_pert_manual %>% filter(organism !="rat") %>% as_tibble()

broad_annot <- read.delim("data/00_db_data/geo_phenotypes_n6470.txt") # from tnatoli
broad2 <- broad_annot %>% 
  group_by(sig_id, class_id) %>%
  dplyr::summarize(list_samples=paste(sample_id, collapse="|"), 
            class_label=paste(unique(class_label), collapse="|"),
            series=paste(unique(series), collapse="|")
            )
broad3 <- broad2 %>% ungroup() %>%
  group_by(sig_id) %>% 
  arrange(class_id) %>%
  dplyr::mutate(class_label=paste(unique(class_label), collapse=";")) %>%
  ungroup() %>%
  pivot_wider(id_cols=c(sig_id, series, class_label), 
              names_from=class_id, values_from=list_samples)
# 6470

creeds_m2 <- creeds_manual %>% 
  select(id, geo_id,  drug_name, drugbank_id, ctrl_ids, pert_ids, organism) # 725
creeds_a2 <- creeds_auto %>% 
  select(id, geo_id,  drug_name, ctrl_ids, pert_ids, organism) # 3817

creeds_a2 <- creeds_a2 %>% 
  mutate(organism=str_replace_all(organism, "\\[|\\]| ", "")) %>%
  separate_rows(organism, sep=",")  %>%
  filter(organism %in% c("human", "mouse"))
  
# CREEDS: applied classifier to n=31905 microarray GEO studies
creeds_a2_gse <- unique(creeds_a2$geo_id) # 308
creeds_m2_gse <- unique(creeds_m2$geo_id) # 299
length(intersect(creeds_a2_gse, creeds_m2_gse)) #6

# Affymetrix studies
broad3_gse <- unique(broad3$series) # 1363
# // TODO: some of these are array express
# many of the are *NOT* drug exposures
length(intersect(creeds_a2_gse, broad3_gse)) # 26
length(intersect(creeds_m2_gse, broad3_gse)) # 61


mapped_gse <- unique((drugbank_mapped %>% filter(study_type=="oligo"))$gse) # 6745 oligo (~8.8k not including)
length(intersect(mapped_gse, creeds_m2_gse)) # 209 (69.9%)
length(intersect(mapped_gse, creeds_a2_gse)) # 134 (43.5%)
length(intersect(mapped_gse, broad3_gse)) # 319 (23.4%)

# why are these different?
#  - different platforms as input data (also we are later)
#  - different performance


# // TODO: map CREEDS auto, BROAD drugs to DrugBank

# <--- alternately, just set up for download CREEDS data ---> #


# write out the mouse/human GSEs

creeds2 <- rbind(creeds_a2 %>% select(geo_id, organism) %>% unique(), 
      creeds_m2 %>% select(geo_id, organism) %>% unique()) %>%
  rename(gse=geo_id) 
creeds2 %>% filter(organism=="human") %>% unique() %>%  # 413
  write_csv("data/01_sample_lists/creeds_human_to_download.csv")

creeds2 %>% filter(organism=="mouse") %>% unique() %>%  # 192
  write_csv("data/01_sample_lists/creeds_mouse_to_download.csv")

# id	cell_type	ctrl_ids	inst_info	GSE	pert_ids	source
# we will look at the different in distributions??


