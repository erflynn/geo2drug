# geo2pubmed_id.R
# E Flynn
# 2/28/2019
#
# Code for mapping GEO IDs to DrugBank IDs:
#  Uses Pubtator + MESH Ids and DrugBank data to accomplish the mapping.
# 
# End result: GSE, DrugBank ID list
#  DrugBank, ATC
#
# Data files:
#    chemical2pubtator
#
# TODO:
#  - update so that this works with PUBMED MESH directly!


require('tidyverse')
require('GEOmetadb')
require('meshr')
require('rjson')

# --- get list of human GSEs, PMIDs ---- #
con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite') 
res <- dbGetQuery(con, "SELECT gse.gse, gpl.gpl, organism, pubmed_id FROM gse JOIN gse_gpl on gse.gse = gse_gpl.gse JOIN gpl ON gse_gpl.gpl=gpl.gpl")
res2 <- separate_rows(res, organism, sep=";\t")
human_data <- filter(res2, organism=="Homo sapiens")
res3 <- filter(human_data, !is.na(pubmed_id)) # 17243 unique GSEs
pubmed_ids <- unique(res3$pubmed_id) # 12816 PMIDs

# --- generate XML commands file --- #
# sink("xml_comamnds.sh")
# sapply(1:13, function(idx){
#   pubmed_st <- paste(pubmed_ids[((idx-1)*1000+1):(idx*1000)], collapse=",")
#   cat(sprintf("efetch -format xml -db pubmed -id \"%s\" > pubmed_%s.xml", pubmed_st, idx))
#   cat("\n")
#   })
# 
# sink()

# create a GSE to pubmed mapping vector
gse_to_pubmed <- res3[!duplicated(res3[,c("gse", "pubmed_id")]),c("gse", "pubmed_id")]
gse_to_pubmed <- rename(gse_to_pubmed, PMID=pubmed_id)

# --- add MeSH ---- #
# read in all the pubtator chemical data, filter by IDs that are associated with a human ID
pubtator <- read.delim("../../drug_expression/drug_labeling/external_data/chemical2pubtator", stringsAsFactors = FALSE)

pubtator_gse <- pubtator[pubtator$PMID %in% pubmed_ids,] # 6654
length(unique(pubtator_gse$MeshID)) # 1800

# fix the MeSH format which is currently MESH:<mesh_id> --> mesh_id
pubtator_gse$MeSH <- sapply(pubtator_gse$MeshID, function(x) 
  {y <- strsplit(x, ":", fixed=TRUE)[[1]]; 
  return(y[[length(y)]])}) # grabs the last in case there are multiple


# link GSE to MESH and write this out
gse_to_mesh <- right_join(gse_to_pubmed, pubtator_gse, by="PMID") # 8655
length(unique(gse_to_mesh$gse)) # 4706 GSEs 
length(unique(gse_to_mesh$Mentions)) # 2468
write.table(gse_to_mesh, file="data/gse_to_mesh.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(data.frame(unique(gse_to_mesh$gse)), file="data/gse_w_drug_mesh.txt", row.names=FALSE, col.names=FALSE, sep="\t")


# condense data where multiple mentions correspond to the same MeSH ID
mesh_to_mention <- select(pubtator_gse, c("MeSH", "Mentions"))
mesh_to_mention <- mesh_to_mention[!duplicated(mesh_to_mention),]
collapsed_mesh <- mesh_to_mention %>% group_by(MeSH) %>% 
  summarise(mention_str=paste(unique(unlist(lapply(Mentions, function(x) tolower(strsplit(x, "\\|")[[1]])))), collapse="|"))
write.table(data.frame(collapsed_mesh$MeSH), file="data/list_mesh.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# note: could collapse further w removing quotes, . , etc

# select the most common mention label per MeSH, this provides an intuitive way to condense!
mesh_labeled_common <- mesh_to_mention %>% group_by(MeSH) %>% 
  summarise(most_common=names(which.max(table(unlist(lapply(Mentions, function(x) tolower(strsplit(x, "\\|")[[1]])))))))

# there are also some CHEBI entries - separate these out
chebi_entries <- sapply(pubtator_gse$MeshID, function(x) strsplit(x, ":")[[1]][[1]]=="CHEBI")
chebi <- pubtator_gse[chebi_entries,] # 275 CHEBI
mesh <-pubtator_gse[!chebi_entries,] # 1525 MeSH


##### ----- START MAPPING ----- #####
 
# map by name
mentions <- separate_rows(collapsed_mesh, mention_str, sep="\\|")
# TODO - double check chemical names are ok

drugbank <- fromJSON(file="data/db_data/drugbank_info.json") 
list.drugs <- lapply(drugbank, function(x) unique(unlist(tolower(c(x$synonyms, x$name)))))
names(list.drugs) <- lapply(drugbank, function(y) y$dbID)
drug.df <- data.frame("drug_names"=unlist(lapply(list.drugs, function(x) paste(x, collapse="|"))))
drug.df$db_id <- names(list.drugs)
drug.df2 <- separate_rows(drug.df, drug_names, sep="\\|")
mapping_tab <- inner_join(mentions, drug.df2, by=c("mention_str"="drug_names"))
dim(mapping_tab) # 576
length(unique(mapping_tab$db_id)) # 486
length(unique(mapping_tab$MeSH)) # 535

list.unmapped <- setdiff(unique(mentions$MeSH), unique(mapping_tab$MeSH)) # 1265

# map by CHEBI
drugbank_chebi <- drugbank[sapply(drugbank, function(x) x$chebi != "")]
chebi_db_df <- data.frame("chebi"=sapply(drugbank_chebi, function(x) x$chebi), 
"db_id"=sapply(drugbank_chebi, function(x) x$dbID), stringsAsFactors = FALSE)
mapping_chebi <- inner_join(chebi, chebi_db_df, c("MeSH"="chebi"))

dim(mapping_chebi) # 136

# - double check CHEBI mapping is the same!
mapped_chebi_and_name <- inner_join(mapping_tab, mapping_chebi, by="MeSH")
mapped_chebi_name_ids <- select(mapped_chebi_and_name, c("MeSH", "db_id.x", "db_id.y"))
mapped_chebi_name_ids <- mapped_chebi_name_ids[!duplicated(mapped_chebi_name_ids),]
table(mapped_chebi_name_ids$db_id.x==mapped_chebi_name_ids$db_id.y)
# 32 match, 2 do not
mapped_chebi_name_ids[mapped_chebi_name_ids$db_id.x!=mapped_chebi_name_ids$db_id.y,]
# 16016 DB03756 DB01775 ### THESE ARE NOT THE SAME
# 36751 DB03088 DB02543 ### THESE ARE SIMILAR BUT NOT THE SAME...
# TODO - resolve these conflicts

# - add drugbank IDs for these to the table
new_chebi <- chebi_db_df[(chebi_db_df$chebi %in% list.unmapped),]  # 22 more mapped
new_chebi <- rename(new_chebi, "MeSH"="chebi")

mapping_tab2 <- rbind(new_chebi, select(mapping_tab, c("MeSH", "db_id")))
mapping_tab2 <- mapping_tab2[!duplicated(mapping_tab2),]
list.unmapped2 <- setdiff(list.unmapped, new_chebi$MeSH)

### - use MESH to map more - ###

# TODO - condense - the MeSH info is in two parts
mesh_info <- fromJSON(file="data/db_data/mesh_info.json") # 233
mesh_data_downloaded <- intersect(mesh$MeSH, names(mesh_info)) # 120
table(mesh_data_downloaded %in% list.unmapped2) # 72 of these are new!

# extract the mappings
mesh_df <- do.call(rbind, 
                   lapply(mesh_info, function(x) 
                     data.frame(lapply(x, function(y) paste(y, collapse=" | ")), stringsAsFactors=FALSE)))
mesh_df$MeSH <- names(mesh_info)
mesh_mapping <- filter(mesh_df, MeSH %in% mesh$MeSH)
# try joining DrugBank on registryNum, CAS 
drugbank_df <- do.call(rbind,
        lapply(drugbank, function(x) data.frame(
  lapply(x[c("synonyms","unii","name","cas","dbID","chebi" )], 
         function(y) paste(y, collapse=" | ")), stringsAsFactors=FALSE)))

mesh_mapping <- rename(mesh_mapping, "cas"="CAS", "unii"="registryNum")
unii_matched <- inner_join(filter(drugbank_df, unii!=""), mesh_mapping, by="unii") # 34
cas_matched <- inner_join(filter(drugbank_df, cas!=""), mesh_mapping, by="cas") # 27

length(intersect(unii_matched$MeSH,cas_matched$MeSH)) # 25 in both
map_unii_cas <- full_join(select(cas_matched, c("dbID", "MeSH")), select(unii_matched, c("dbID", "MeSH")), by="MeSH")
filter(map_unii_cas, dbID.x!=dbID.y) # ok so none conflict - yay!

cas_unii <- rbind(select(cas_matched, c("dbID", "MeSH")), select(unii_matched, c("dbID", "MeSH")))
cas_unii <- cas_unii[!duplicated(cas_unii),]

head(mapping_tab2)
cas_unii <- rename(cas_unii, db_id=dbID)
table(cas_unii$MeSH %in% list.unmapped2) # only 3 are unmapped..., 33 already
mapping_tab3 <- rbind(mapping_tab2, cas_unii[,colnames(mapping_tab2)])
mapping_tab3 <- mapping_tab3[!duplicated(mapping_tab3),] # 563
write.table(mapping_tab3, file="data/mesh_db_mapping_0228.txt", row.names=FALSE, sep="\t", quote=FALSE)



still.unmapped <- setdiff(list.unmapped2, names(mesh_info)) # 1171 <-- download THESE!
write.table(data.frame(still.unmapped), file="data/mesh_to_download.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# --- load the second mesh file ---#
mesh_info2 <- fromJSON(file="data/db_data/mesh_info2.json")
length(mesh_info2)
mesh_df2 <- do.call(rbind, 
                   lapply(mesh_info2, function(x) 
                     data.frame(lapply(x, function(y) paste(y, collapse=" | ")), stringsAsFactors=FALSE)))
mesh_df2$MeSH <- names(mesh_info2)
head(mesh_df2)
table(mesh_df2$name == "Add") # 207 of these....
mesh.failed <- filter(mesh_df2, name=="Add")$MeSH
mesh.failed2 <- mesh.failed[sapply(mesh.failed, function(x) substr(x, 1,1) %in% c("C", "D"))] # only 19 are actual MeSH Ids - these do not have anything on the webpages

mesh_present2 <- filter(mesh_df2, name!="Add")


mesh_mapping2 <- filter(mesh_present2, MeSH %in% mesh$MeSH)
mesh_mapping2 <- rename(mesh_mapping2, "cas"="CAS", "unii"="registryNum")
unii_matched2 <- inner_join(filter(drugbank_df, unii!=""), mesh_mapping2, by="unii") # 55
cas_matched2 <- inner_join(filter(drugbank_df, cas!=""), mesh_mapping2, by="cas") # 24

length(intersect(unii_matched2$MeSH,cas_matched2$MeSH)) # 22 in both
map_unii_cas2 <- full_join(select(cas_matched2, c("dbID", "MeSH")), select(unii_matched2, c("dbID", "MeSH")), by="MeSH")
filter(map_unii_cas2, dbID.x!=dbID.y) # ok so none conflict - yay!


cas_unii2 <- rbind(select(cas_matched2, c("dbID", "MeSH")), select(unii_matched2, c("dbID", "MeSH")))
cas_unii2 <- cas_unii2[!duplicated(cas_unii2),]

cas_unii2 <- rename(cas_unii2, db_id=dbID)
table(cas_unii2$MeSH %in% list.unmapped2) # 57!
mapping_tab4 <- rbind(mapping_tab3, cas_unii2[,colnames(mapping_tab3)])
mapping_tab4 <- mapping_tab4[!duplicated(mapping_tab4),] # 563

## ---- map synonyms ---- ##
mesh_df_full <- rbind(mesh_df, mesh_present2) # 964

syn_map <- mesh_df_full[,c("synonyms", "name", "MeSH")]
syn_map2 <- unite(syn_map, "names", synonyms, name, sep="|")

syn_map3 <- separate_rows(syn_map2, names, sep="\\|")
drug_syn <- drugbank_df[,c("synonyms", "name", "dbID")]
drug_syn2 <- unite(drug_syn, "names", synonyms, name, sep="|")
drug_syn3 <- separate_rows(drug_syn2, names, sep="\\|")

# clean up the strings - remove spaces, tolower
drug_syn3$names <- sapply(drug_syn3$names, function(x) trimws(tolower(x)))
syn_map3$names <- sapply(syn_map3$names, function(x) trimws(tolower(x)))
mesh_to_db_syn <- inner_join(drug_syn3, syn_map3, by="names")
mesh_to_db_syn <- mesh_to_db_syn[!duplicated(mesh_to_db_syn),] # 208 rows, 136 MeSH IDs

mesh_syn_new <- select(mesh_to_db_syn[!mesh_to_db_syn$MeSH %in% mapping_tab4$MeSH,] , c("MeSH", "dbID"))
mesh_syn_new <- mesh_syn_new[!duplicated(mesh_syn_new),] # 47
mesh_syn_new <- rename(mesh_syn_new, "db_id"="dbID")
mapping_tab5 <- rbind(mapping_tab4, mesh_syn_new[colnames(mapping_tab4)])


list.unmapped.final <- setdiff(unique(mentions$MeSH), unique(mapping_tab5$MeSH)) #1149, 961 with MeSH IDs
mesh.unmapped.final <- list.unmapped.final[sapply(list.unmapped.final, function(x) substr(x, 1,1) %in% c("C", "D"))]

write.table(data.frame(list.unmapped.final), file="data/unmapped_mesh.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# what do we do with all these that don't map to drugbank??
View(mesh_labeled_common[mesh_labeled_common$MeSH %in%list.unmapped.final,])
View(mesh_df_full[mesh_df_full$MeSH %in% list.unmapped.final,])

mesh_unmapped <- mesh_df_full[mesh_df_full$MeSH %in% list.unmapped.final,]
head(mesh_unmapped) # 942


# map these 6... then give up for now
new_mapping <- mesh_unmapped[(mesh_unmapped$registryNum %in% drugbank_df$cas),][,c("MeSH", "registryNum")]
new_mapping <- rename(new_mapping, "cas"="registryNum")
new_by_cas <- inner_join(new_mapping, drugbank_df[,c("cas", "dbID")],by="cas")
new_by_cas <- rename(new_by_cas, "db_id"="dbID")
mapping_tab6 <- rbind(mapping_tab5, new_by_cas[colnames(mapping_tab4)])
write.table(mapping_tab6, file="data/mesh_db_mapping_0302.txt", row.names=FALSE, sep="\t", quote=FALSE) # 673 map
gse_mesh_db <- full_join(gse_mesh, mapping_tab6, by="MeSH")
write.table(gse_mesh_db, file="data/gse_mesh_db.txt", row.names=FALSE, sep="\t", quote=FALSE)
