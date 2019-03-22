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

creeds_data <- rbind(select(drug_pert_manual, colnames(drug_pert_auto)), drug_pert_auto)
creeds_gse <- unique(creeds_data$geo_id) # 343

# ---- load the BROAD annotations ---- #
broad_annot <- read.delim2("data/ref_data/geo_phenotypes_n6470.txt") # 193,370 samples
broad_gse <- unique(broad_annot$series) # 1363

# reformat the BROAD annotations s.t. the GSMs for the perts/ctls follow a similar format to the CREEDS data
pert_annot <- broad_annot %>% group_by(sig_id) %>% filter(class_id=="A") %>% summarise("GSE"=unique(series), pert_ids = paste(sample_id, collapse="|"), pert_label=paste(unique(class_label), collapse="|") )
ctl_annot <- broad_annot %>% group_by(sig_id) %>% filter(class_id=="B") %>% summarise(ctl_ids = paste(sample_id, collapse="|"), ctl_label=paste(unique(class_label), collapse="|") )
broad_comb_annot <- full_join(pert_annot, ctl_annot, by="sig_id")
broad_comb_annot$source <- "broad_manual"
dim(broad_comb_annot) # 6470 trt/ctl instances

table(sapply(broad_comb_annot$GSE, function(x) str_detect(x, "GSE"))) # 58 are arrayExpress, 1305 are GSE

# TODO - put together CREEDS and BROAD data
# BROAD data doesn't have a "drug_name" field


# ---- load the annotations we put together ---- #
gse_mesh_db <- read.delim("data/gse_mesh_db.txt")

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

# TODO - mesh to drugbank mapping
# TODO - cell line info

# NOW LOOK AT THE OUTPUT

