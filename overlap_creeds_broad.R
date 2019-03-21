# overlap_creeds_broad.R 
# E Flynn
# 3/5/2019
#
# Look at the overlap with CREEDS, BROAD annotations
# Then put this together:
#   GSE, ctl GSMs, pert GSMs, DrugBank ID, DrugName, sex labels?, tissue/cell line labels?

require('tidyverse')

# load the CREEDS data
drug_pert_auto <- read.csv("../../drug_expression/drug_labeling/data/single_drug_perturbations-p1.0.csv", stringsAsFactors = FALSE)
colnames(drug_pert_auto)
drug_pert_auto$source <- "creeds_auto"

drug_pert_manual <- read.csv("../../drug_expression/drug_labeling/data/single_drug_perturbations-v1.0.csv", stringsAsFactors = FALSE)
drug_pert_manual$source <- "creeds_manual"
rbind(select(drug_pert_manual, colnames(drug_pert)))
creeds_gse <- unique(drug_pert_manual$geo_id, drug_pert_auto$geo_id) # 343
# load the BROAD annotations
broad_annot <- read.delim2("data/geo_phenotypes_n6470.txt")
dim(broad_annot)
length(unique(broad_annot$class_label))

broad_gse <- unique(broad_annot$series) # 1363

# reformat this into a collapsed table of trt vs contol
head(broad_annot)




# count number of controls, treatments, sex label?
length(intersect(broad_gse,creeds_gse)) # 61



# load the annotations we put together
mapping_tab6 <-read.table("data/mesh_db_mapping_0302.txt", header=TRUE)
gse_mesh <- read.delim("data/gse_to_mesh.txt", header=TRUE)
gse_mesh_db <- full_join(gse_mesh, mapping_tab6, by="MeSH")
head(gse_mesh_db)
mesh_gse <- unique(gse_mesh_db$gse) # 4707

overlapping_b <- intersect(mesh_gse, broad_gse) # 359
overlapping_c <- intersect(creeds_gse, mesh_gse) # 77

# this is very little overlap - why are we getting so little?


# look at an overlapping case 
gse_likely_trt_ctl <- gse_mesh_db[gse_mesh_db$gse %in% union(overlapping_b, overlapping_c),]
unique(gse_likely_trt_ctl$db_id)

pert_annot <- broad_annot %>% group_by(sig_id) %>% filter(class_id=="A") %>% summarise("GSE"=unique(series), pert_ids = paste(sample_id, collapse="|"), pert_label=paste(unique(class_label), collapse="|") )
ctl_annot <- broad_annot %>% group_by(sig_id) %>% filter(class_id=="B") %>% summarise(ctl_ids = paste(sample_id, collapse="|"), ctl_label=paste(unique(class_label), collapse="|") )
comb_annot <- full_join(pert_annot, ctl_annot, by="sig_id")
dim(comb_annot) # 6470 groups

table(sapply(comb_annot$GSE, function(x) str_detect(x, "GSE"))) # 58 are arrayExpress  1305 are GSE
head(comb_annot)

comb_annot_filt <- filter(comb_annot, GSE %in% overlapping_b)
broad_w_db <- inner_join(comb_annot, select(gse_mesh_db, "gse", "Mentions", "db_id"), by=c("GSE"="gse"))
head(broad_w_db)

# grab drug, tissue, sex info

label_mat <- read.csv("data/label_mat.csv")
head(label_mat)
label_mat <- rename(label_mat, "gsm" ="X")
ale_data <- read.csv("../../drug_expression/drug_labeling/ale_processing/ale_combined_data.csv")
comb_labels <- left_join(label_mat, select(ale_data, c("gsm", "gse", "gpl", "text_sex", "text_tissue_name", "cell_line")))          




# NOW LOOK AT THE OUTPUT



