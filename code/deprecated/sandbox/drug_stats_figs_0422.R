
options(stringsAsFactors=FALSE)
require('tidyverse')
cl <- read.delim("data/labeled_data/cell_line_mapped_gse.txt", sep=" ")
drug <- read.delim("data/labeled_data/drugbank_mapped_gse.txt", sep=" ")

drug2 <- filter(drug, ! name %in% c("Glucose", "D-glucose", "Oxygen", "Nitrogen", "L-Glutamine", "Sucrose"))
length(unique(drug2$dbID)) # [1] 1087
length(unique(drug2$gse)) # [1] 8480
drug_counts <- table(drug2$name)
head(drug_counts[order(-drug_counts)], 60)

drug_cl <- left_join(drug2, cl, by="gse")

# add in organism data!
require('GEOmetadb')
formatted_list_gses <- paste(drug_cl$gse, collapse="\', \'")
con <- dbConnect(SQLite(),'../../drug_expression/dataset_identification/GEOmetadb.sqlite') # TODO move
study_mapping <- dbGetQuery(con, 
                            sprintf("SELECT gse, gpl.gpl, gpl.organism FROM gse_gpl JOIN gpl ON gpl.gpl = gse_gpl.gpl WHERE gse IN ('%s');", formatted_list_gses))
dbDisconnect(con)

# UGH SEPARATE THE STUPID ROWS

# FILTER for mouse/rat only
org_dat <- filter(study_mapping, organism %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus"))
drug_plus_org <- inner_join(drug_cl, org_dat)

head(drug_plus_org)
drug_plus_org2 <- drug_plus_org %>% mutate(is_cell_line=(cell_line | !is.na(accession)))
 # metadata sex
table(drug_plus_org2[,c("organism", "is_cell_line")])

drug_counts <- table(drug2$name)
head(drug_counts[order(-drug_counts)], 60)

mouse_counts <- read.csv("../sex_labeling/data/mouse_rat/mouse_gse_counts.csv") %>% mutate("organism"="mouse")
rat_counts <- read.csv("../sex_labeling/data/mouse_rat/rat_gse_counts.csv")  %>% mutate("organism"="rat")
human_counts <- read.csv("../sex_labeling/data/mouse_rat/human_gse_counts.csv")  %>% mutate("organism"="human")
mrh_counts <- do.call(rbind, list(mouse_counts, rat_counts, human_counts))
study_w_counts <- left_join(drug_plus_org2, mrh_counts, by="gse")
dim(study_w_counts)
# cell line sex
cell_info_df <- read.csv("data/db_data/cellosaurus_df.txt")
cell_sex <- cell_info_df[,c("sex", "accession")]
cell_sex$accession <- sapply(cell_sex$accession, tolower)
comb_dat <- left_join(study_w_counts, cell_sex, by="accession")

table(study_dedup[,c("organism.y", "study_type")])
table(study_dedup$organism.y)
study_dedup <- study_w_counts[,c("study_type", "organism.y", "is_cell_line", "gse")] %>% unique()
table(study_dedup[,c("organism.y", "is_cell_line")])
table(study_dedup$is_cell_line)


## <---- PIE CHARTS, breakdown by organism of labels ----> ##

# look at CREEDS/BROAD


# what have I actually labeled?

# what does the smoking data look like?