# compareDrugNames.R
# E Flynn
# 

require('tidyverse')

# --- 0. MeSH vs Pubtator --- #
pmid_to_mesh <- fromJSON(file="data/pmid_to_mesh.json") # 11771, 1491
pmid <- as.character(gse_to_pubmed[gse_to_pubmed$gse=="GSE68895",]$PMID)
pmid_to_mesh[pmid]
filter(pubtator_gse, PMID==pmid)
# so it appears that the pubtator annotations have a lot more info than the PMID to MeSH
# --> check this out in more detail: what is the overlap?
head(length(names(pmid_to_mesh)))
pubtator_to_mesh <- pubtator_gse[,c("PMID", "MeSH")]
# remove the CHEBI??
pubtator_pmid_to_mesh <- split(pubtator_to_mesh$MeSH, pubtator_to_mesh$PMID)
length(pubtator_pmid_to_mesh) # 3651

# look at the overlap of the ones that overlap??
length(intersect(names(pmid_to_mesh), names(pubtator_pmid_to_mesh))) # 3470 - most of them?



# --- 1. Compare MeSH/Pubtator mentions with drugbank --- #

gse_mesh_db <- read.delim("data/gse_mesh_db.txt")
head(gse_mesh_db)



# --- 1.5. Use the CREEDS manual data as a sort of gold standard --- #

drug_manual <- drug_pert_manual[,c("drug_name", "drugbank_id", "geo_id")]
drug_manual <- rename(drug_manual, "gse"="geo_id") %>% unique()
manual2 <- left_join(drug_manual, unique(drug[,c("gse", "dbID", "name")]))

manual3 <- manual2 %>% group_by(gse) %>% summarize(true_name=paste(unique(drug_name), collapse="|"), 
                                        true_id=paste(unique(drugbank_id), collapse="|"), 
                                        pred_id=paste(unique(dbID), collapse="|"), 
                                        pred_name=paste(unique(name), collapse="|"),
                                        overlapping_ids=paste(intersect(dbID, drugbank_id), collapse="|"),
                                        num_overlap=length(intersect(dbID, drugbank_id)),
                                        num_true=length(unique(drugbank_id)), 
                                        unlabeled=(pred_id=="NA"))
manual3 <- is.na(manual3$pred_id)
table(manual3[,c("unlabeled", "num_overlap")])
# 170/390


# --- 2. Label the BROAD/CREEDS data --- #



comb_c_b <- read.csv("data/broad_creeds_annot_combined.csv")
head(comb_c_b)
cl <- read.delim("data/labeled_data/cell_line_mapped_gse.txt", sep=" ")
drug <- read.delim("data/labeled_data/drugbank_mapped_gse.txt", sep=" ")

cell <- left_join(comb_c_b, cl, by=c("GSE"="gse"))
labeled_c_b <- left_join(cell, drug, by=c("GSE"="gse"))
head(labeled_c_b)

labeled_mesh <- left_join(comb_c_b, gse_mesh_db, by=c("GSE"="gse"))

all_annot_hc <- inner_join(labeled_c_b, gse_mesh_db, by=c("GSE"="gse"))
both_annot <- filter(all_annot_hc, !is.na(db_id) & !is.na(dbID)) %>% select("GSE", "inst_info", "name", "dbID", "db_id") %>% unique()
overlap <- both_annot %>% group_by(GSE) %>% summarize(db_overlap=paste(intersect(dbID, db_id), collapse="|"))
View(filter(both_annot, GSE %in% filter(overlap,db_overlap=="")$GSE))