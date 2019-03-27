# geo2pubmed_id.R
# E Flynn
# 2/28/2019
#
# Code for mapping GEO IDs to MESH IDs:
#  Uses Pubtator + MESH Ids and DrugBank data to accomplish the mapping.
# 
# End result: GSE, MESH list
#
# Data files from other directories:
#    chemical2pubtator
#    GEOMetaDB

require('tidyverse')
require('GEOmetadb')
require('rjson')
source("drug_count_utils.R")

# --- get list of human GSEs, PMIDs ---- #
con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite') 
res <- dbGetQuery(con, "SELECT gse.gse, gpl.gpl, organism, pubmed_id FROM gse JOIN gse_gpl on gse.gse = gse_gpl.gse JOIN gpl ON gse_gpl.gpl=gpl.gpl")
res2 <- separate_rows(res, organism, sep=";\t")
human_data <- filter(res2, organism=="Homo sapiens")
res3 <- filter(human_data, !is.na(pubmed_id)) # 17243 unique GSEs
pubmed_ids <- unique(res3$pubmed_id) # 12816 PMIDs

# create a GSE to pubmed mapping vector
gse_to_pubmed <- res3[!duplicated(res3[,c("gse", "pubmed_id")]),c("gse", "pubmed_id")]
gse_to_pubmed <- rename(gse_to_pubmed, PMID=pubmed_id)

# --- generate XML commands file --- #
# sink("xml_comamnds.sh")
# sapply(1:13, function(idx){
#   pubmed_st <- paste(pubmed_ids[((idx-1)*1000+1):(idx*1000)], collapse=",")
#   cat(sprintf("efetch -format xml -db pubmed -id \"%s\" > pubmed_%s.xml", pubmed_st, idx))
#   cat("\n")
#   })
# 
# sink()


# --- add MeSH ---- #
# MeSH mapping was done by downloading the PUBMED XML and Pubtator annotations for the PMIDs 
#  This was then filtered for drug-related MeSH IDs!
#
# load the mesh data (after 'parse_pubmed_xml.py' has been run)
pmid_to_mesh <- fromJSON(file="data/db_data/pmid_to_mesh.json")
pmid_mesh_df <- vec_to_df(pmid_to_mesh)
colnames(pmid_mesh_df) <- c("MeSH", "PMID")
pmid_mesh_df$PMID <- sapply(pmid_mesh_df$PMID, as.numeric)
length(unique(inner_join(gse_to_pubmed, pmid_mesh_df)$gse))

pmid_mesh_df_long <- separate_rows(pmid_mesh_df, MeSH, sep="\\|")
all_mesh_ids <- unique(pmid_mesh_df_long$MeSH) # 1491
mesh_type <- sapply(all_mesh_ids, grabMeshType)
table(mesh_type) 
# D    Q 
# 1424   67 
# <--- why are there are only D, Q - what about other categories??? 
pmid_mesh_df_long2 <- filter(pmid_mesh_df_long, grabMeshType(MeSH)=="D")

# --- read in all the pubtator chemical data, filter by IDs that are associated with a human ID --- #
pubtator <- read.delim("../../drug_expression/drug_labeling/external_data/chemical2pubtator", stringsAsFactors = FALSE)

pubtator_gse <- pubtator[pubtator$PMID %in% pubmed_ids,] # 6654
length(unique(pubtator_gse$MeshID)) # 1800

length(unique(inner_join(gse_to_pubmed, pubtator_gse)$gse))

# fix the MeSH format which is currently MESH:<mesh_id> --> mesh_id
pubtator_gse$MeSH <- sapply(pubtator_gse$MeshID, function(x) 
  {y <- strsplit(x, ":", fixed=TRUE)[[1]]; 
  return(y[[length(y)]])}) # grabs the last in case there are multiple

# condense data where multiple mentions correspond to the same MeSH ID
mesh_to_mention <- pubtator_gse[, c("MeSH", "Mentions")]
mesh_to_mention <- mesh_to_mention[!duplicated(mesh_to_mention),]
collapsed_mesh <- mesh_to_mention %>% group_by(MeSH) %>% 
  summarise(mention_str=paste(unique(unlist(lapply(Mentions, function(x) tolower(strsplit(x, "\\|")[[1]])))), collapse="|"))
#write.table(data.frame(collapsed_mesh$MeSH), file="data/list_mesh.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# note: could collapse further w removing quotes, . , etc

# select the most common mention label per MeSH, this provides an intuitive way to condense!
mesh_labeled_common <- mesh_to_mention %>% group_by(MeSH) %>% 
  summarise(most_common=names(which.max(table(unlist(lapply(Mentions, function(x) tolower(strsplit(x, "\\|")[[1]])))))))

# there are also some CHEBI entries - separate these out
chebi_entries <- sapply(pubtator_gse$MeshID, function(x) strsplit(x, ":")[[1]][[1]]=="CHEBI")
pubtator_chebi <- pubtator_gse[chebi_entries,] # 275 CHEBI
pubtator_mesh <- pubtator_gse[!chebi_entries,] # 1525 MeSH
pubtator_mesh2 <- pubtator_mesh[,c("PMID", "MeSH")]
pubtator_mesh2 <- pubtator_mesh2[!duplicated(pubtator_mesh2),]

# ---- put together pubmed XML and pubtator data ---- #
pubtator_gse2 <- pubtator_gse[,c("MeSH", "PMID")]
pubtator_gse2 <- pubtator_gse2[!duplicated(pubtator_gse2),]
comb_mesh_annot <- rbind(pubtator_gse2, pmid_mesh_df_long2[!duplicated(pmid_mesh_df_long2),])
table(duplicated(comb_mesh_annot)) # only 182 pairs overlap... remaining 18,361 do not! this is interesting...

all_mesh_ids2 <- union(unique(pubtator_mesh$MeSH), all_mesh_ids[mesh_type=="D"]) # 2853

# OUTPUT:
#  PMID GSE MESHID
comb_mesh_annot$PMID <- sapply(comb_mesh_annot$PMID, as.numeric)
gse_mesh <- inner_join(gse_to_pubmed, comb_mesh_annot)
length(unique(gse_mesh$MeSH))
table(sapply(unique(gse_mesh$MeSH), grabMeshType) %in% c("D", "C")) 
length(unique(gse_mesh$gse))
write.table(gse_mesh, file="data/gse_to_mesh.txt", row.names=FALSE, sep="\t")