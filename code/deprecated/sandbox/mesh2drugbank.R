# mesh2drugbank.R
# E Flynn
# 3/27/2019
#
# Code for mapping MESH IDs to Drugbank
#
# 1. Map by Name/Synonym
#    from MeSH name/synyonyms <--> DrugBank name/synonyms
#


# ---- 1. Map by name/synonym ---- #
# load vocabs
drug_vocab <- read.delim("data/db_data/drugbank_vocab.txt")
mesh_vocab <- read.delim("data/db_data/mesh_vocab_w_db.txt")

# remove things that are less than 3 characters - these are gonna be too ambiguous
mesh_vocab2 <- mesh_vocab[sapply(mesh_vocab$name, nchar) >= 3, ]
drug_vocab2 <- drug_vocab[sapply(drug_vocab$name, nchar) >= 3, ]

# match the names
mesh_drugbank <- inner_join(mesh_vocab2, drug_vocab2) 

length(unique(mesh_drugbank$MeSH)) # 509
length(unique(mesh_drugbank$dbID)) # 519

# load the mentions and use these to map
collapsed_mesh <- read.delim("data/db_data/pubtator_mesh_to_mention.txt")
mesh_mentions <- separate_rows(collapsed_mesh, mention_str, sep="\\|")
mesh_mentions <- rename(mesh_mentions, "name"="mention_str")
mesh_mentions2 <- mesh_mentions[sapply(mesh_mentions$name, nchar) >= 3, ]
mentions_drugbank <- inner_join(mesh_mentions2, drug_vocab2) # 

# combine and check for multi-mapping
combined_mapping <- rbind(mesh_drugbank, mentions_drugbank)
combined_mapping <- combined_mapping[!duplicated(combined_mapping),] # 656 MeSH to 603 drugbank, 1034 records

combined2 <- combined_mapping[!duplicated(combined_mapping[,c("MeSH", "dbID")]),]




# ---- 2. Map by CHEBI/UNII/CAS IDs ---- #


# load drugbank data
drugbank_df <- read.delim("data/db_data/drugbank_parsed.txt")

# load MeSH data
mesh_df_full <- read.delim("data/db_data/mesh_info_df.txt")


# map by CHEBI
drugbank_chebi <- drugbank_df[!is.na(drugbank_df$chebi),]
drugbank_chebi$chebi <- sapply(drugbank_chebi$chebi, as.character)
mapping_chebi <- inner_join(collapsed_mesh, drugbank_chebi, c("MeSH"="chebi"))

dim(mapping_chebi) # 56 out of 275

# try unii, CAS 
mesh_mapping <- rename(mesh_df_full, "cas"="CAS", "unii"="registryNum")
unii_matched <- inner_join(filter(drugbank_df, unii!=""), mesh_mapping, by="unii") # 281
cas_matched <- inner_join(filter(drugbank_df, cas!=""), mesh_mapping, by="cas") # 190

length(intersect(unii_matched$MeSH,cas_matched$MeSH)) # 174 in both
map_unii_cas <- full_join(cas_matched[,c("dbID", "MeSH")], unii_matched[,c("dbID", "MeSH")], by="MeSH")
filter(map_unii_cas, dbID.x!=dbID.y) # one conflicts :/ 

cas_unii <- rbind(cas_matched[,c("dbID", "MeSH")], unii_matched[,c("dbID", "MeSH")])
cas_unii <- cas_unii[!duplicated(cas_unii),]



# map a few with mixed up unii/CAS
mixed_up <- mesh_df_full[mesh_df_full$registryNum!= "",]
new_mapping <- rename(mixed_up[,c("MeSH", "registryNum")], "cas"="registryNum")
new_by_cas <- inner_join(new_mapping, drugbank_df,by="cas")

mixed_up2 <- mesh_df_full[mesh_df_full$CAS!= "",]
new_mapping2 <- rename(mixed_up2[,c("MeSH", "CAS")], "unii"="CAS")
new_by_unii <- inner_join(new_mapping2, drugbank_df,by="unii")


# ---- 3. Put together and look for conflicts ---- #
all_mapping <- do.call(rbind, list(
        combined_mapping[,c("MeSH", "dbID")],
  mapping_chebi[,c("MeSH", "dbID")],
  cas_unii[,c("MeSH", "dbID")],
  new_by_cas[,c("MeSH", "dbID")],
  new_by_unii[,c("MeSH", "dbID")]))

all_mapping2 <- all_mapping[!duplicated(all_mapping),] # 702 MeSH to 637 DrugBank
length(unique(all_mapping2$MeSH))
length(unique(all_mapping2$dbID))
mesh_to_db <- split(all_mapping2$dbID, all_mapping2$MeSH)
multi_map_mesh <- mesh_to_db[sapply(mesh_to_db, length) > 1] # 21
db_to_mesh <- split(all_mapping2$MeSH, all_mapping2$dbID)
multi_map_db <- db_to_mesh[sapply(db_to_mesh, length) > 1] # 68

# TODO - disambiguate multi-mapping 


write.table(all_mapping2, "data/mesh_to_drugbank_0327.txt", row.names=FALSE, sep="\t")

# Which still are not mapped?
# <-- issue many of these are "Pharmacology"
drug_likely_mesh <- unique(drug_likely$MESHID)

missing <- setdiff(drug_likely_mesh,all_mapping2$MeSH)
View(unique(drug_likely[drug_likely$MESHID %in% missing,]$MESHTERM))
head(collapsed_mentions)

mesh_all <- rbind(mesh_mentions2[,colnames(mesh_vocab2)], mesh_vocab2)
mesh_all2 <- filter(mesh_all, MeSH %in% missing)
#head(drug_vocab2[drug_vocab2$dbID=="DB00783",]) # estrogen doesn't map b/c "estradiol"
#mesh_all2[mesh_all2$name=="estrogen",]
#mesh_all[mesh_all$MeSH=="D004967",]
