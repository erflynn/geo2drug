# mesh2drugbank.R
# E Flynn
# 3/27/2019
#
# Code for mapping MESH IDs to Drugbank


require('MeSH.db')


gse_to_mesh <- read.delim("data/gse_to_mesh.txt")
all_chem_ids <- unique(gse_to_mesh$MeSH)
all_mesh_ids <- all_chem_ids[sapply(all_chem_ids, grabMeshType) %in% c("D", "C")] 

mesh_db_info <- select(MeSH.db, keys=all_mesh_ids, columns=columns(MeSH.db), keytype="MESHID") # SLOW
# also this ONLY grabs 1912 out of 2853 - where are the other 941
mesh_not_in_db <- setdiff(all_mesh_ids, mesh_db_info$MESHID)

mesh_df_full <- read.delim2("data/db_data/mesh_info_df.txt")
length(intersect(mesh_db_info$MESHID, mesh_df_full$MeSH))

write.table(data.frame(mesh_no_info), file="data/mesh_to_download_v3.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)


######### WHAT DATA CORRESPONDS TO DRUGS???? ###########
# likely a drug if it contains a "pharm" ID
drug_likely <- filter(mesh_db_info, QUALIFIER=="pharmacology") 
drug_likely <- drug_likely[,c("MESHID", "MESHTERM", "SYNONYM")]
drug_likely_mesh_id <- unique(drug_likely$MESHID) # 1122
drug_unlikely <- setdiff(mesh_db_info$MESHID, drug_likely_mesh_id) # 790

# --- try to map to drugbank based on name/synonym --- #

# put together names + synonyms / mentions for MESH data
mesh_db_data <- mesh_db_info %>% mutate(SYNONYM=sapply(SYNONYM, function(x) strsplit(x, "\\|")[[1]][[1]])) # parse junk from synonym
synonym_data <- mesh_db_data[,c("MESHID", "SYNONYM")]
name_data <- mesh_db_data[,c("MESHID", "MESHTERM")]
name_data <- name_data[!duplicated(name_data),] # condense to function
name_data <- rename(name_data, "SYNONYM"="MESHTERM")
name_syn_df <- rbind(name_data, synonym_data)
name_syn_df$SYNONYM <- sapply(name_syn_df$SYNONYM, tolower)

mentions <- separate_rows(collapsed_mesh, mention_str, sep="\\|")
mentions2 <- rename(mentions, "drug_names"="mention_str")
name_syn_df2 <- rename(name_syn_df, "MeSH"="MESHID", "drug_names"="SYNONYM")
mesh_syn_comb <- rbind(mentions2, name_syn_df2)

# load drugbank data
drugbank_df <- read.delim("data/db_data/drugbank_parsed.txt")


# map by name/synonym
mapping_tab <- inner_join(mesh_syn_comb, drug.df2)
length(unique(mapping_tab$db_id)) # 552
length(unique(mapping_tab$MeSH)) # 597
length(setdiff(mapping_tab$MeSH, drug_likely_mesh_id)) # 256 
length(union(mapping_tab$MeSH, drug_likely_mesh_id)) # 1378
list.unmapped <- setdiff(all_mesh_ids2, unique(mapping_tab$MeSH)) # 2321

# map by CHEBI
drugbank_chebi <- drugbank[sapply(drugbank, function(x) x$chebi != "")]
chebi_db_df <- data.frame("chebi"=sapply(drugbank_chebi, function(x) x$chebi), 
                          "db_id"=sapply(drugbank_chebi, function(x) x$dbID), stringsAsFactors = FALSE)
mapping_chebi <- inner_join(pubtator_chebi, chebi_db_df, c("MeSH"="chebi"))

dim(mapping_chebi) # 137

# - double check CHEBI mapping is the same!
mapped_chebi_and_name <- inner_join(mapping_tab, mapping_chebi, by="MeSH")
mapped_chebi_name_ids <- mapped_chebi_and_name[, c("MeSH", "db_id.x", "db_id.y")]
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

mapping_tab2 <- rbind(new_chebi, mapping_tab[,c("MeSH", "db_id")])
mapping_tab2 <- mapping_tab2[!duplicated(mapping_tab2),]
list.unmapped2 <- setdiff(list.unmapped, new_chebi$MeSH)


### - use MESH to map more - ###

# try joining DrugBank on registryNum, CAS 
mesh_mapping <- rename(mesh_mapping, "cas"="CAS", "unii"="registryNum")
unii_matched <- inner_join(filter(drugbank_df, unii!=""), mesh_mapping, by="unii") # 89
cas_matched <- inner_join(filter(drugbank_df, cas!=""), mesh_mapping, by="cas") # 51

length(intersect(unii_matched$MeSH,cas_matched$MeSH)) # 174 in both
map_unii_cas <- full_join(cas_matched[,c("dbID", "MeSH")], unii_matched[,c("dbID", "MeSH")], by="MeSH")
filter(map_unii_cas, dbID.x!=dbID.y) # one conflicts :/ 

cas_unii <- rbind(cas_matched[,c("dbID", "MeSH")], unii_matched[,c("dbID", "MeSH")])
cas_unii <- cas_unii[!duplicated(cas_unii),]


cas_unii <- rename(cas_unii, db_id=dbID)
table(cas_unii$MeSH %in% list.unmapped2) # 41

mapping_tab3 <- rbind(mapping_tab2, cas_unii[,colnames(mapping_tab2)])
mapping_tab3 <- mapping_tab3[!duplicated(mapping_tab3),] # 676


## ---- map synonyms ---- ##

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

mesh_syn_new <-mesh_to_db_syn[!mesh_to_db_syn$MeSH %in% mapping_tab3$MeSH,]
mesh_syn_new <- mesh_syn_new[, c("MeSH", "dbID")]
mesh_syn_new <- mesh_syn_new[!duplicated(mesh_syn_new),] # 47
mesh_syn_new <- rename(mesh_syn_new, "db_id"="dbID")
mapping_tab5 <- rbind(mapping_tab3, mesh_syn_new[colnames(mapping_tab3)])


list.unmapped.final <- setdiff(list.unmapped, unique(mapping_tab5$MeSH)) #2260, 976 with MeSH IDs
mesh.unmapped.final <- list.unmapped.final[sapply(list.unmapped.final, function(x) substr(x, 1,1) %in% c("C", "D"))]
mesh.unmapped.final_d <- setdiff(mesh.unmapped.final, drug_unlikely) # 1470
unmapped.not.downloaded <- setdiff(mesh.unmapped.final_d, mesh_df_full$MeSH) # 457
to.download <- union(unmapped.not.downloaded, mesh_missing) # 642 --> 320
write.table(data.frame(mesh_missing), "data/mesh_to_download_v2.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(data.frame(list.unmapped.final), file="data/unmapped_mesh_0322.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# what do we do with all these that don't map to drugbank??
#View(mesh_labeled_common[mesh_labeled_common$MeSH %in%list.unmapped.final,])
#View(mesh_df_full[mesh_df_full$MeSH %in% list.unmapped.final,])

mesh_unmapped <- mesh_df_full[mesh_df_full$MeSH %in% mesh.unmapped.final_d,] # 1013

# map these 6... then give up for now
new_mapping <- mesh_unmapped[(mesh_unmapped$unii %in% drugbank_df$cas & mesh_unmapped$unii != ""),][,c("MeSH", "unii")]
new_mapping <- rename(new_mapping, "cas"="unii")
new_by_cas <- inner_join(new_mapping, drugbank_df[,c("cas", "dbID")],by="cas")
new_by_cas <- rename(new_by_cas, "db_id"="dbID")

mapping_tab6 <- rbind(mapping_tab5, new_by_cas[colnames(mapping_tab5)]) # 686
write.table(mapping_tab6, file="data/mesh_db_mapping_0322.txt", row.names=FALSE, sep="\t", quote=FALSE) # 704 map
# -- hmmm I was expecting more... -- #


# filter
gse_mesh_db <- full_join(gse_mesh, mapping_tab6, by="MeSH")
write.table(gse_mesh_db, file="data/gse_mesh_db.txt", row.names=FALSE, sep="\t", quote=FALSE)





# what are the GSEs associated with these data
gse_to_pubmed$PMID <- sapply(gse_to_pubmed$PMID, as.character)
gse_mapping_d <- inner_join(gse_to_pubmed, filter(comb_mesh_annot, MeSH %in% drug_likely_mesh_id))
head(gse_mapping_d)
length(unique(gse_mapping_d$gse)) # 6,549
length(unique(gse_mapping_d$PMID)) # 4,933
# RESULT:
# <GSE> : <PMID> : <MESH>
#   for MESH that are likely drug

# ok wait --> actually filtering again...





# link GSE to MESH and write this out
gse_to_mesh <- right_join(gse_to_pubmed, pubtator_gse, by="PMID") # 8655
length(unique(gse_to_mesh$gse)) # 4706 GSEs 
length(unique(gse_to_mesh$Mentions)) # 2468
write.table(gse_to_mesh, file="data/gse_to_mesh.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(data.frame(unique(gse_to_mesh$gse)), file="data/gse_w_drug_mesh.txt", row.names=FALSE, col.names=FALSE, sep="\t")





# ----- MAP TO DRUGBANK ----- #





