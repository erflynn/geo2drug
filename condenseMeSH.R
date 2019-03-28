# condenseMeSH.R
# E Flynn
# 2/27/2019
# 
# Code for condensing MeSH data into a single dataframe

require('rjson')
require('tidyverse')
require('MeSH.db')
source('drug_count_utils.R')

mesh_info <- fromJSON(file="data/db_data/mesh_info.json")
mesh_info2 <- fromJSON(file="data/db_data/mesh_info2.json")
mesh_info3 <- fromJSON(file="data/db_data/mesh_info3.json")
mesh_info4 <- fromJSON(file="data/db_data/mesh_info4.json")

# extract the mappings
mesh_df <- vec_to_df(mesh_info) %>% rename("MeSH"="id")
mesh_df2 <- vec_to_df(mesh_info2) %>% rename("MeSH"="id")
mesh_df3 <- vec_to_df(mesh_info3) %>% rename("MeSH"="id")
mesh_df4 <- vec_to_df(mesh_info4) %>% rename("MeSH"="id")


mesh_df_full <- do.call(rbind, list(mesh_df, mesh_df2, mesh_df3, mesh_df4))

table(mesh_df_full$name == "Add") # 275 of these....
mesh.failed <- filter(mesh_df_full, name=="Add")$MeSH
mesh.failed2 <- mesh.failed[sapply(mesh.failed, grabMeshType) %in% c("C", "D")] # only 19 are actual MeSH Ids - these do not have anything on the webpages

mesh_df_full <- filter(mesh_df_full, name!="Add")
write.table(mesh_df_full, "data/db_data/mesh_info_df.txt", sep="\t", row.names=FALSE)

# create the vocabulary
mesh_name_syn <- createNameSynVocab(mesh_df_full, "MeSH")
write.table(mesh_name_syn, file="data/db_data/mesh_vocab.txt", sep="\t", row.names=FALSE)

# ------- other MESH data ------ #
gse_to_mesh <- read.delim("data/gse_to_mesh.txt")
all_chem_ids <- unique(gse_to_mesh$MeSH)
all_mesh_ids <- all_chem_ids[sapply(all_chem_ids, grabMeshType) %in% c("D", "C")] 

mesh_db_info <- select(MeSH.db, keys=all_mesh_ids, columns=columns(MeSH.db), keytype="MESHID") # SLOW
# also this ONLY grabs 1912 out of 2853 - where are the other 941
mesh_not_in_db <- setdiff(all_mesh_ids, mesh_db_info$MESHID)
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
mesh_db_data <- rename(mesh_db_data, "synonyms"="SYNONYM", "name"="MESHTERM", "MeSH"="MESHID")
mesh_db_syn <- createNameSynVocab(mesh_db_data, "MeSH")
mesh_vocab <- rbind(mesh_db_syn, mesh_name_syn)
mesh_vocab <- mesh_vocab[!duplicated(mesh_vocab),]

write.table(mesh_vocab, file="data/db_data/mesh_vocab_w_db.txt", row.names=FALSE, sep="\t") ## 22781

