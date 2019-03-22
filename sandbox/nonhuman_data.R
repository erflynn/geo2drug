# nonhuman_data.R
# E Flynn
# 3/5


res <- dbGetQuery(con, "SELECT gse.gse, gpl.gpl, organism, pubmed_id FROM gse JOIN gse_gpl on gse.gse = gse_gpl.gse JOIN gpl ON gse_gpl.gpl=gpl.gpl")
res2 <- separate_rows(res, organism, sep=";\t")
non_human <- filter(res2, organism!="Homo sapiens" & !is.na(pubmed_id)) # 31,771
length(unique(non_human$pubmed_id)) # 20,056
length(unique(non_human$gse)) # 26,978

non_h_pmid_df <- filter(non_human, pubmed_id %in% nonh_drug_pmids)
length(unique(non_h_pmid_df$gse)) # 8159 GSEs

org_breakdown <- table(non_h_pmid_df$organism)
head(org_breakdown[order(-org_breakdown)], 30)
org_breakdown <- data.frame(org_breakdown[order(-org_breakdown)])
head(org_breakdown, 15)
# what about the non-human data?
#  GSE to MeSH to drugbank for these


pubtator <- read.delim("../../drug_expression/drug_labeling/external_data/chemical2pubtator", stringsAsFactors = FALSE)
nonh_drug_pmids <- intersect(non_human$pubmed_id, pubtator$PMID)
length(nonh_drug_pmids) # 6488 PMIDs with drug data
pubtator_nonh <- filter(pubtator,PMID %in% nonh_drug_pmids)
length(unique(pubtator_nonh$MeshID)) # 3022 mesh IDs

# start mapping to drugbank
pubtator_nonh$MeSH <- sapply(pubtator_nonh$MeshID, function(x) 
{y <- strsplit(x, ":", fixed=TRUE)[[1]]; 
return(y[[length(y)]])})
head(pubtator_nonh)
length(intersect(pubtator_nonh$MeSH, mapping_tab6$MeSH )) # at least 428 have DrugBank IDs
# TODO:
# - finish map to drugbank

head(non_h_pmid_df)
non_h_pmid_df2 <- rename(non_h_pmid_df, "PMID"="pubmed_id")
nonh_mapped <- inner_join(non_h_pmid_df2, pubtator_nonh)
nonh_mapped2 <- left_join(nonh_mapped, mapping_tab6)