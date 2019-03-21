# geo2pubmed_id.R
# E Flynn
# 2/28/2019
#
# Code for mapping GEO IDs to DrugBank IDs:
#  Uses Pubtator + MESH Ids to accomplish the mapping.
#
# End result: GSE, DrugBank ID list
#  DrugBank, ATC

require('tidyverse')
require('GEOmetadb')
require('meshr')
require('rjson')

con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite') 

# TaxonID = 9606
res <- dbGetQuery(con, "SELECT gse.gse, gpl.gpl, organism, pubmed_id FROM gse JOIN gse_gpl on gse.gse = gse_gpl.gse JOIN gpl ON gse_gpl.gpl=gpl.gpl")
res2 <- separate_rows(res, organism, sep=";\t")
human_data <- filter(res2, organism=="Homo sapiens")
res3 <- filter(human_data, !is.na(pubmed_id))
length(unique(res3$gse)) # 17243
pubmed_ids <- unique(res3$pubmed_id) # 12816

gse_to_pubmed <- res3[!duplicated(res3[,c("gse", "pubmed_id")]),c("gse", "pubmed_id")]
gse_to_pubmed <- rename(gse_to_pubmed, PMID=pubmed_id)

# read in all the pubtator chemical data
pubtator <- read.delim("../../drug_expression/drug_labeling/external_data/chemical2pubtator", stringsAsFactors = FALSE)
head(pubtator)

pubtator_gse <- pubtator[pubtator$PMID %in% pubmed_ids,] # 6654
length(unique(pubtator_gse$MeshID)) # 1800
pubtator_gse$MeSH <- sapply(pubtator_gse$MeshID, function(x) 
  {y <- strsplit(x, ":", fixed=TRUE)[[1]]; 
  return(y[[length(y)]])})


# condense the mapping to mentions
mesh_to_mention <- select(pubtator_gse, c("MeSH", "Mentions"))
mesh_to_mention <- mesh_to_mention[!duplicated(mesh_to_mention),]

collapsed_mesh <- mesh_to_mention %>% group_by(MeSH) %>% 
  summarise(mention_str=paste(unique(unlist(lapply(Mentions, function(x) tolower(strsplit(x, "\\|")[[1]])))), collapse="|"))

# could collapse further w removing quotes, . , etc
write.table(data.frame(collapsed_mesh$MeSH), file="data/list_mesh.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# alternate <-- select the most common
mesh_labeled_common <- mesh_to_mention %>% group_by(MeSH) %>% summarise(most_common=names(which.max(table(unlist(lapply(Mentions, function(x) tolower(strsplit(x, "\\|")[[1]])))))))

# filter chebi vs mesh entries
chebi_entries <- sapply(pubtator_gse$MeshID, function(x) strsplit(x, ":")[[1]][[1]]=="CHEBI")
chebi <- pubtator_gse[chebi_entries,] #  275 unique
mesh <-pubtator_gse[!chebi_entries,] #1525 unique

length(unique(pubtator_gse$PMID)) # 3587

gse_to_mesh <- right_join(gse_to_pubmed, pubtator_gse, by="PMID") # 8655
length(unique(gse_to_mesh$gse)) # 4706 GSEs - should I download ALL of these? 
length(unique(gse_to_mesh$Mentions)) # 2468
write.table(gse_to_mesh, file="data/gse_to_mesh.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(data.frame(unique(gse_to_mesh$gse)), file="data/gse_w_drug_mesh.txt", row.names=FALSE, col.names=FALSE, sep="\t")
##### ----- START MAPPING ----- #####

# Steps:
#   MeSH, DrugBank ID (or other)
#
# Eventually:
#   GSE, PMID, DrugBank IDs



# MAP BY NAME 

mentions <- separate_rows(collapsed_mesh, mention_str, sep="\\|")
# TODO - double check chemical names are ok

drugbank <- fromJSON(file="data/drugbank_info.json") 
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

mesh_info <- fromJSON(file="data/mesh_info.json") # 233
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




# TODO:
# check already existing
# - check chebi mapping

# map the new ones

still.unmapped <- setdiff(list.unmapped2, names(mesh_info)) # 1171 <-- download THESE!
write.table(data.frame(still.unmapped), file="data/mesh_to_download.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# --- load the second mesh file ---#
mesh_info2 <- fromJSON(file="data/mesh_info2.json")
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
mesh_unmapped
mesh_unmapped[mesh_unmapped$name=="Immunoferon",]

# map these 6... then give up lol
new_mapping <- mesh_unmapped[(mesh_unmapped$registryNum %in% drugbank_df$cas),][,c("MeSH", "registryNum")]
new_mapping <- rename(new_mapping, "cas"="registryNum")
new_by_cas <- inner_join(new_mapping, drugbank_df[,c("cas", "dbID")],by="cas")
new_by_cas <- rename(new_by_cas, "db_id"="dbID")
mapping_tab6 <- rbind(mapping_tab5, new_by_cas[colnames(mapping_tab4)])
write.table(mapping_tab6, file="data/mesh_db_mapping_0302.txt", row.names=FALSE, sep="\t", quote=FALSE) # 673 map
head(mapping_tab6)

##### ----------------------- ####
##### ----------------------- ####
##### ----------------------- ####

##### ----------------------- OLD
# 
# 
# ### ----  write out the commands for download with eutilities ---- ###
# sink("xml_comamnds.sh")
# sapply(1:13, function(idx){
#   pubmed_st <- paste(pubmed_ids[((idx-1)*1000+1):(idx*1000)], collapse=",")
#   cat(sprintf("efetch -format xml -db pubmed -id \"%s\" > pubmed_%s.xml", pubmed_st, idx))
#   cat("\n")
#   })
# 
# sink()
# 
# ### ---- end ---- ###
# 
# ### ---- use python to parse the XML ---- ###
# # python parse_pubmed_xml.py
# 
# ### --- load the data --- ###
# 
# require('rjson')
# 
# pmid_to_mesh <- fromJSON(file="data/pmid_to_mesh.json") # 11771, 1491
# mesh_ids <- unique(unlist(pmid_to_mesh))
# chemical_mesh <- mesh_ids[sapply(mesh_ids, function(x) substr(x, 1, 1)=="D")] # 1494
# 
# 
# 
# pubtator$MeSH <- sapply(pubtator$MeshID, function(x) {y <- strsplit(as.character(x), ":", fixed=TRUE)[[1]]; ifelse(length(y)==2, y[[2]], y)})
# 
# filt_pubtator <- pubtator[pubtator$MeSH %in% chemical_mesh,]
# list_mentions <- unique(filt_pubtator$Mentions)
# 
# # a lot of these don't appear to be unique - and many are multiple names sorted
# length(unique(filt_pubtator$MeshID)) # only 233 of these
# 
# # --> COLLAPSE!
# 
# head(collapsed_mesh)
# dim(collapsed_mesh)
# 
# # query for MESH
# 
# 
# # map these to drugbank/ATC
# drugbank <- fromJSON(file="drugbank_info.json") 
# 
# 
# 
# list.drugs <- lapply(drugbank, function(x) unique(unlist(tolower(c(x$synonyms, x$name)))))
# names(list.drugs) <- lapply(drugbank, function(y) y$dbID)
# 
# 
# 
# intersect(unique(mesh_labeled_common$most_common), unlist(list.drugs)) # only 50 intersect out of 233 :(
# setdiff(unique(mesh_labeled_common$most_common), unlist(list.drugs))
# 
# # try to map
# 
# head(res3)
# head(pmid_to_mesh)
# head(filt_pubtator) # pubtator filtered by the PMID to MESH list
# # not sure this was the correct way to go...
# # why did I do this additional mapping
# # TODO - double check that these match!!
# filt_pubtator$PMID <- sapply(filt_pubtator$PMID, as.character)
# length(intersect(names(pmid_to_mesh), filt_pubtator$PMID)) # 754
# length(names(pmid_to_mesh)) # 11,771
# length(unique(filt_pubtator$PMID)) # 1,143,649
# 
# overlapping <- intersect(names(pmid_to_mesh), filt_pubtator$PMID)
# missing <- setdiff(names(pmid_to_mesh), filt_pubtator$PMID)
# 
# 
# 
# # now filter the PUBMED mapping
# # which have MESH IDs?
# 
# # MESH to ATCCC??? (or use UMLS)
# 
# 
# # alternately - try to figure out how to map using names/synonyms to drugbank
# #  I think this could actually work?
# 
# # write out a list of GEO to PMID to MESH ID
# #  - start sex labeling all of these!
# 
# 
# 
# 
# 
# # write out a list of PMID to ATCCC
# 
