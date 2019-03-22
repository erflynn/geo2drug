# refine_drug_annot.R
# E Flynn
# 

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

# look at the set diff??
 

# what would be better validation?
#   CREEDS?? BROAD data?
