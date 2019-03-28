# overlap_creeds_broad.R 
# E Flynn
# 3/5/2019
#
# Look at the overlap with CREEDS, BROAD annotations
# Then put this together:
#   GSE, ctl GSMs, pert GSMs, DrugBank ID, DrugName, sex labels?, tissue/cell line labels?

require('tidyverse')
options(stringsAsFactors=FALSE)
# ----- load the CREEDS data ----- #
drug_pert_auto <- read.csv("../../drug_expression/drug_labeling/data/single_drug_perturbations-p1.0.csv", stringsAsFactors = FALSE)
drug_pert_auto$source <- "creeds_auto"

drug_pert_manual <- read.csv("../../drug_expression/drug_labeling/data/single_drug_perturbations-v1.0.csv", stringsAsFactors = FALSE)
drug_pert_manual$source <- "creeds_manual"
# note - the manual data has additional fields including the DrugBank ID
# TODO - check overlap w drugbank IDs

creeds_data <- rbind(drug_pert_manual[, colnames(drug_pert_auto)], drug_pert_auto)
creeds_gse <- unique(creeds_data$geo_id) # 709

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

table(sapply(broad_comb_annot$GSE, function(x) str_detect(x, "GSE"))) # 355 are arrayExpress, 6115 are GSE

# --- put together CREEDS and BROAD data --- #
# BROAD data doesn't have a "drug_name" field
head(broad_comb_annot)
head(creeds_data)
broad_comb_annot <- rename(broad_comb_annot, "inst_info"="pert_label", "id"="sig_id", "ctrl_ids"="ctl_ids")
broad_comb_annot$cell_type <- NA
creeds_data <- rename(creeds_data, "inst_info"="drug_name", "GSE"="geo_id")
overlapping_cols <- intersect(colnames(creeds_data), colnames(broad_comb_annot))
comb_c_b <- rbind(creeds_data[,overlapping_cols], broad_comb_annot[,overlapping_cols])
write.csv(comb_c_b, file="data/broad_creeds_annot_combined.csv", row.names=FALSE)

hc_drug_gse <- union(broad_gse, creeds_gse) # 1987

# --- map these to DRUGBANK --- #

drugbank_data <- read.delim2("data/db_data/drugbank_parsed.txt")

drug_name_syn <- read.delim2("data/db_data/drugbank_vocab.txt")

# broad data
broad_comb_annot$inst_info <-sapply(broad_comb_annot$inst_info, tolower)

# --- divide into drug exposures vs. other exposures vs. disease! --- #
trt_inst <- str_detect(broad_comb_annot$inst_info, "treat|agent|before|after|expose")
oe_inst <- str_detect(broad_comb_annot$inst_info, "overexp|over-exp")
kd_inst <- str_detect(broad_comb_annot$inst_info, "knock|kd|deplete")
normal_inst <- str_detect(broad_comb_annot$inst_info, "normal|control|healthy") # a couple of these are treated
rna_inst <- str_detect(broad_comb_annot$inst_info, "transfect|mir|sirna|transduce")
#View(broad_comb_annot[normal_inst, c("inst_info", "ctl_label")])
#View(broad_comb_annot[!(trt_inst | oe_inst | kd_inst | normal_inst | rna_inst), c("inst_info", "ctl_label")])

#require('fuzzyjoin')

# from https://stackoverflow.com/questions/32914357/dplyr-inner-join-with-a-partial-string-match
partial_join <- function(x, y, by_x, pattern_y){
  idx_x <- sapply(y[[pattern_y]], grep, fixed(x[[by_x]]))
idx_y <- sapply(seq_along(idx_x), function(i) rep(i, length(idx_x[[i]])))

df <- dplyr::bind_cols(x[unlist(idx_x), , drop = F],
                       y[unlist(idx_y), , drop = F])
return(df)
}


# need to escape all the characters in a chemical mapping 
escapeALLChars <- function(x){
  tmp <- gsub( "\\{", "\\\\{", x)
  tmp <- gsub( "\\}", "\\\\}", tmp)
  tmp <- gsub( "\\[", "\\\\[", tmp)
  tmp <- gsub( "\\]", "\\\\]", tmp)
  tmp <- gsub( "\\(", "\\\\(", tmp)
  tmp <- gsub( "\\)", "\\\\)", tmp)
  return(tmp)
}

drug_name_syn$name <- sapply(drug_name_syn$name, escapeALLChars)
# filter out names that are smaller than 5 characters - these will match too many things
drug_name_syn2 <- drug_name_syn[(sapply(drug_name_syn$name, function(x) nchar(x) >= 5)),]
# go back and match unmatched w/ 3-5 char??


broad_mapped <- partial_join(broad_comb_annot, drug_name_syn2, by_x="inst_info", pattern_y="name")

# filter by super common
head(table(broad_mapped$name)[order(-table(broad_mapped$name))], 50)
broad_mapped <- filter(broad_mapped, !name %in% c("helium", "cyclo", "carbon")) # common and non-specific

max_mapping <- broad_mapped %>% group_by(inst_info) %>% filter(name==name[which.max(nchar(name))]) # extract the lONGEST MATCH
max_mapping <- max_mapping[!duplicated(max_mapping),] # 1293 out of 6470 - this is not great...
head(table(max_mapping$name)[order(-table(max_mapping$name))], 40)
# filter by TF/IDF?

# what about the ones that don't map?
# -- first: filter out RNAs, etc
not_rna <- broad_comb_annot[!( oe_inst | kd_inst | normal_inst | rna_inst),] # 3966
head(setdiff(not_rna$inst_info, max_mapping$inst_info), 20)
#View(broad_comb_annot[!broad_comb_annot$id %in% max_mapping$id,c("inst_info")])
# do any of them map with a 3/4-letter mapping?
mapping2 <- partial_join(filter(broad_comb_annot, !id %in% broad_mapped$id), drug_name_syn, by_x="inst_info", pattern_y="name")
head(table(broad_mapped$name)[order(-table(broad_mapped$name))], 50)



# manual have drugbank IDs

# automatic do not - these labels are effectively CRAP
 # 0. <-- do these map thru MESH?
 # count how many have an assoc PMID + MESH ID
 # 1. download metadata and try to extract???
drug_pert_auto2 <- separate_rows(drug_pert_auto, drug_name, sep="\\|")
auto_mapped <- inner_join(drug_pert_auto2, drug_vocab2, c("drug_name"="name"))
dim(auto_mapped)
length(unique(auto_mapped$id)) # 2496 -- this is actually a fair number!

drug_pert_manual2 <- rename(drug_pert_manual, "dbID"="drugbank_id")
creeds_data2 <- rbind(drug_pert_manual2[, colnames(auto_mapped)], auto_mapped)


# BROAD look interesting - I think we need to run str_detect on all of these?
head(broad_comb_annot$inst_info)
broad_mapped <- rename(broad_mapped, "geo_id"="GSE") 
broad_mapped <- rename(broad_mapped, "drug_name"="name")
overlapping_cols <- intersect(colnames(broad_mapped), colnames(creeds_data2))
combined_df_drugbank <- rbind(creeds_data2[,overlapping_cols], broad_mapped[,overlapping_cols])
dim(combined_df_drugbank)

write.table(combined_df_drugbank, file="data/hc_combined_df_drugbank.txt", sep="\t", row.names=FALSE)

# which are human
require('GEOmetadb')
con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite') 
res <- dbGetQuery(con, "SELECT gse.gse, gpl.gpl, organism, pubmed_id FROM gse JOIN gse_gpl on gse.gse = gse_gpl.gse JOIN gpl ON gse_gpl.gpl=gpl.gpl")
res2 <- separate_rows(res, organism, sep=";\t") 
human_data <- filter(res2, organism=="Homo sapiens")
organism_labels <- select(res2, "gse", "organism")
hc_drug_plus_org <- filter(organism_labels, gse %in% hc_drug_gse)
write.csv(hc_drug_plus_org, file="data/hc_drug_org.csv", row.names=FALSE)
human_hc_drug <- intersect(hc_drug_gse, human_data$gse) # 1608

arrayexp <- setdiff(hc_drug_gse, hc_drug_plus_org$gse) # 58 ARRAY EXPRESS studies
write.csv(data.frame(arrayexp), file="data/hc_array_exp.csv", row.names=FALSE)
# ---- load the annotations we put together ---- #
gse_mesh_db <- read.delim("data/gse_mesh_db.txt")

gse_to_download <- setdiff(human_hc_drug, gse_mesh_db$gse) # 1149
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
# TODO - cell line info

# ---- create a "count per drug exposure" table ---- #


# this is like the count.per.study table BUT has counts divided by pert_id, ctl_id
# num_f_p, num_f_c, num_m_p, num_m_c, study_type



# ---- remake plots looking at these data ---- #

# barplot of sex breakdown of GSEs


# ATC breakdown of GSEs



