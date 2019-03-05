# examine_drug_freq.R
# E Flynn
# 3/4/2019
#
# Code for examining drug frequency in GEO.
#   - what are the ATC codes?
#   - plot the common drugs

require('tidyverse')
require('rjson')
mapping_tab6 <-read.table("data/mesh_db_mapping_0302.txt", header=TRUE)

gse_mesh <- read.delim("data/gse_to_mesh.txt", header=TRUE)

gse_mesh_db <- full_join(gse_mesh, mapping_tab6, by="MeSH")
head(gse_mesh_db)

# add drugbank names + ATC codes
drugbank <- fromJSON(file="drugbank_info.json") 
drugbank_df <- do.call(rbind,
                       lapply(drugbank, function(x) data.frame(
                         lapply(x[c("synonyms","unii","name","cas","dbID","chebi" , "ATC")], 
                                function(y) paste(y, collapse=" | ")), stringsAsFactors=FALSE)))


# what are the most frequent drugs?
name_df <- drugbank_df[,c("dbID", "name", "ATC")]
name_df <- rename(name_df, "db_id"="dbID")
gse_mesh_db_n <- left_join(gse_mesh_db, name_df)
gse_mesh_db2 <- gse_mesh_db_n[!is.na(gse_mesh_db_n$db_id),]
gse_mesh_db_na <- gse_mesh_db_n[is.na(gse_mesh_db_n$db_id),]
drug.counts <- table(as.factor(gse_mesh_db2$name))
head(drug.counts[order(-drug.counts)], 30)
require('ggplot2')
gse_mesh_db2$name <- factor(gse_mesh_db2$name, names(drug.counts)[order(-drug.counts)])
ggplot(gse_mesh_db2, aes(name)) + geom_histogram(stat="count", bins=30)+
  ylab("Number of GSEs")+ggtitle("Distribution of drug frequencies in GEO Studies")+
  xlab("Drugbank Drugs (n=600)")+theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())
#
# look at ones with drugbank ids
# hmm appears that many are estrogen, androgen, etc
# TODO
#  go through and map these
summary(as.factor(gse_mesh_db_na$Mentions))

# breakdown by drug class
atc_classes <- sapply(gse_mesh_db2$ATC, function(x) ifelse(x=="", NA, substr(x, 1, 1)))
table(is.na(atc_classes)) # 1489
atc_counts <- data.frame(table(atc_classes[!is.na(atc_classes)]))
colnames(atc_counts) <- c("class", "count")
atc_names <- read.delim("data/atc_classes.txt", header=FALSE)
head(atc_names)
colnames(atc_names) <- c("class", "descript")
atc_info <- inner_join(atc_counts, atc_names)
head(atc_info)
atc_info$descript <- factor(atc_info$descript, atc_info$descript)
ggplot(atc_info, aes(x=class, y=count, fill=descript))+geom_histogram(stat="identity")+ggtitle("ATC breakdown")

gse_w_study_info <- left_join(gse_mesh_db2, count.per.study)
head(gse_w_study_info)
gse_w_study_info$class <- sapply(gse_w_study_info$ATC, function(x) ifelse(x=="", NA, substr(x, 1, 1)))
gse_w_study_info2 <- left_join(gse_w_study_info, atc_names)
ggplot(filter(gse_w_study_info2, !is.na(study_type) & !is.na(descript)), aes(x=class, fill=descript))+geom_histogram(stat="count")+ggtitle("ATC breakdown")+facet_grid(study_type ~ .)+theme(text = element_text(size=20))
head(data.frame("count"=summary(as.factor(filter(gse_w_study_info2, study_type=="f")$name))), 10)
head(data.frame("count"=summary(as.factor(filter(gse_w_study_info2, study_type=="m")$name))), 10)
head(data.frame("count"=summary(as.factor(filter(gse_w_study_info2, study_type=="both")$name))), 10)

atc_tiss <- inner_join(gse_w_study_info2, study_tiss_count)
atc_tiss <- left_join(atc_tiss, atc_info[,c("class", "descript")])
atc_tiss2 <-filter( atc_tiss, class != "L")
tiss_class <- table(atc_tiss2[,c("descript", "common_tiss")])
tiss_class/rowSums(tiss_class)
require('gplots')
heatmap.2(tiss_class, col=colorRampPalette(c("white", "blue", "darkblue"))(n=500),  
          Rowv=FALSE, Colv=FALSE, trace="none", 
          margins=c(6,22))
filter(atc_tiss, class=="V" & common_tiss=="Muscle" & !is.na(gse))
#  GSE35764 fructose load
#  GSE9105 insulin
table(inner_join(gse_w_study_info2[gse_w_study_info2$name=="Estradiol",], study_tiss_count)[,c("study_type", "common_tiss")])
table(inner_join(gse_w_study_info2[gse_w_study_info2$name=="Fluorouracil",], study_tiss_count)[,c("study_type", "common_tiss")])

filter(atc_tiss, name=="Fluorouracil")
select(filter(atc_tiss, name=="Doxorubicin"), c("gse", "PMID", "MeshID", "name", "num_f", "num_m", "study_type", "common_tiss"))
table(filter(atc_tiss, name=="Doxorubicin")[,c("study_type", "common_tiss")])


count_table <- table(gse_w_study_info[,c("study_type", "class")])
chisq.p <- lapply(colnames(count_table), function(x)
  chisq.test(count_table[c("f", "m"),x])$p.value)
names(chisq.p) <- colnames(count_table)
chisq.p[chisq.p < 0.05/ncol(count_table)]

require('GEOmetadb')
con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite') 
res <- dbGetQuery(con, "SELECT gse.gse, gpl.gpl, organism, pubmed_id, gse.submission_date FROM gse JOIN gse_gpl on gse.gse = gse_gpl.gse JOIN gpl ON gse_gpl.gpl=gpl.gpl")
res2 <- filter(res, organism=="Homo sapiens")
res3 <- filter(res2, !is.na(pubmed_id))
length(unique(res3$gse)) # 16904
pubmed_ids <- unique(res3$pubmed_id) # 12576

# look at date patterns
drug_w_date <- left_join(select(gse_w_study_info2, c("gse", "PMID", "name", "class", "study_type")), select(res3, c("gse", "submission_date")))
drug_w_date$Month <-as.Date(cut(as.Date(drug_w_date$submission_date), breaks = "month"))
drug_w_date$Year <-as.Date(cut(as.Date(drug_w_date$submission_date), breaks = "year"))
year_counts <- filter(drug_w_date, !is.na(study_type)) %>% group_by(Year, study_type) %>% summarize(num_per_year=n())

ggplot(year_counts, aes(x=Year, y=num_per_year, fill=study_type)) + geom_histogram(stat="identity", position="dodge")+ylab("Number of studies per year")+theme(text = element_text(size=20))


# what about the non-human data?
#  GSE to MeSH to drugbank for these
non_human <- filter(res,organism!="Homo sapiens" & !is.na(pubmed_id)) # 31,771
length(unique(non_human$pubmed_id)) # 20,056
length(unique(non_human$gse)) # 26,978

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


# counts
non_h_pmid_df <- filter(non_human, pubmed_id %in% nonh_drug_pmids)
length(unique(non_h_pmid_df$gse)) # 8159 GSEs

org_breakdown <- table(non_h_pmid_df$organism)
head(org_breakdown[order(-org_breakdown)], 30)
org_breakdown <- data.frame(org_breakdown[order(-org_breakdown)])
head(org_breakdown, 15)
# TODO 
# - some studies have human + other - need to make sure we are getting these!


# TODO
# - how does this relate to year?
# ideas:
#  (are there any drugs with interesting trends?)
head(non_h_pmid_df)
non_h_pmid_df2 <- rename(non_h_pmid_df, "PMID"="pubmed_id")
nonh_mapped <- inner_join(non_h_pmid_df2, pubtator_nonh)
nonh_mapped2 <- left_join(nonh_mapped, mapping_tab6)


# TODO
# - how does it overlap with other data?
#  - manual data from CREEDS, the BROAD






