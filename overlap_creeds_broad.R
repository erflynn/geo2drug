# overlap_creeds_broad.R 
# E Flynn
# 3/5/2019
#
# Look at the overlap with CREEDS, BROAD annotations
# Then put this together:
#   GSE, ctl GSMs, pert GSMs, DrugBank ID, DrugName, sex labels?, tissue/cell line labels?

require('tidyverse')

# ----- load the CREEDS data ----- #
drug_pert_auto <- read.csv("../../drug_expression/drug_labeling/data/single_drug_perturbations-p1.0.csv", stringsAsFactors = FALSE)
drug_pert_auto$source <- "creeds_auto"

drug_pert_manual <- read.csv("../../drug_expression/drug_labeling/data/single_drug_perturbations-v1.0.csv", stringsAsFactors = FALSE)
drug_pert_manual$source <- "creeds_manual"
# note - the manual data has additional fields including the DrugBank ID
# TODO - check overlap w drugbank IDs

creeds_data <- rbind(drug_pert_manual[, colnames(drug_pert_auto)], drug_pert_auto)
creeds_gse <- unique(creeds_data$geo_id) # 343

# ---- load the BROAD annotations ---- #
broad_annot <- read.delim2("data/ref_data/geo_phenotypes_n6470.txt") # 193,370 samples
broad_gse <- unique(broad_annot$series) # 1363

# reformat the BROAD annotations s.t. the GSMs for the perts/ctls follow a similar format to the CREEDS data
pert_annot <- broad_annot %>% group_by(sig_id) %>% filter(class_id=="A") %>% 
  summarise("GSE"=unique(series), pert_ids = paste(sample_id, collapse="|"), pert_label=paste(unique(class_label), collapse="|") )
ctl_annot <- broad_annot %>% group_by(sig_id) %>% filter(class_id=="B") %>% 
  summarise(ctl_ids = paste(sample_id, collapse="|"), ctl_label=paste(unique(class_label), collapse="|") )
broad_comb_annot <- full_join(pert_annot, ctl_annot, by="sig_id")
broad_comb_annot$source <- "broad_manual"
dim(broad_comb_annot) # 6470 trt/ctl instances

table(sapply(broad_comb_annot$GSE, function(x) str_detect(x, "GSE"))) # 58 are arrayExpress, 1305 are GSE

# TODO - put together CREEDS and BROAD data
# BROAD data doesn't have a "drug_name" field

hc_drug_gse <- union(broad_gse, creeds_gse) # 1987

# which are human
require('GEOmetadb')
con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite') 
res <- dbGetQuery(con, "SELECT gse.gse, gpl.gpl, organism, pubmed_id FROM gse JOIN gse_gpl on gse.gse = gse_gpl.gse JOIN gpl ON gse_gpl.gpl=gpl.gpl")
res2 <- separate_rows(res, organism, sep=";\t")
human_data <- filter(res2, organism=="Homo sapiens")

human_hc_drug <- intersect(hc_drug_gse, human_data$gse) # 1608

# ---- load the annotations we put together ---- #
gse_mesh_db <- read.delim("data/gse_mesh_db.txt")


gse_to_download <- setdiff(hc_drug_gse, gse_mesh_db$gse)
write.table(data.frame(gse_to_download), file="data/gse_to_dowload_creeds_broad.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collapse by study --> each study contains multiple mentions and their drugbank IDs
#  doing this because we will be joining on the study
gse_mesh_collapse <- select(gse_mesh_db, "gse", "Mentions", "db_id")
gse_mesh_collapse <- gse_mesh_collapse %>% group_by(gse) %>% summarise("Mentions"=paste(Mentions, collapse="|"), "db_ids"=paste(db_id, collapse="|")) 
mesh_gse <- gse_mesh_collapse$gse # 4707

# what is the overlap at the GSE level
overlap_c_b <- intersect(broad_gse,creeds_gse) # 61
overlapping_b <- intersect(mesh_gse, broad_gse) # 359
overlapping_c <- intersect(creeds_gse, mesh_gse) # 77

# Hmmm... this is very little overlap

# TODO - problem - we only downloaded the overlapping data --> need to download the rest
#  creeds_gse
#

# look at annotations that overlap
broad_w_db <- inner_join(broad_comb_annot, gse_mesh_collapse, by=c("GSE"="gse")) # 1572 drug trt/ctl inst
head(broad_w_db[,c("pert_label", "Mentions")], 20) # some of these look reasonable, some look like BS

creeds_w_db <- inner_join(creeds_data, gse_mesh_collapse, by=c("geo_id"="gse")) # 611 drug trt/ctl inst
head(creeds_w_db[,c("drug_name", "Mentions")], 20) # again some look like BS but others do check out


# ---- grab drug, tissue, sex, cell line info ---- #
label_mat <- read.csv("data/expr_label_mat.csv") # expression-based labels from three methods + tissue
label_mat <- rename(label_mat, "gsm" ="X")
ale_data <- read.csv("../../drug_expression/drug_labeling/ale_processing/ale_combined_data.csv") # text-based labels
comb_labels <- left_join(label_mat, select(ale_data, c("gsm", "gse", "gpl", "text_sex", "text_tissue_name", "cell_line")))          

# filter the comb_labels by what is in creeds, broad

drugbank_df <- read.delim("data/db_data/drugbank_parsed.txt") # drugbank info

# TODO - cell line info

# ---- create a "count per drug exposure" table ---- #


# this is like the count.per.study table BUT has counts divided by pert_id, ctl_id
# num_f_p, num_f_c, num_m_p, num_m_c, study_type



# ---- remake plots looking at these data ---- #

# barplot of sex breakdown of GSEs


# ATC breakdown of GSEs



