# searchMetadataDrugNames.R
# E Flynn
# 4/14/2019

require('tidyverse')
require('tidytext')
options(stringsAsFactors = FALSE)

# read in the metadata
# ---- GSE data ---- #
gse_data <- read.csv("data/db_data/gse_all_geo_info.csv")
gse_data$str <-sapply(1:nrow(gse_data), function(i) {paste(gse_data[i,c("title", "summary", "overall_design")][!is.na(gse_data[i,c("title", "summary", "overall_design")])], collapse=" ")}) 
#gsm_filt_w_gse <- fread("data/db_data/gsm_all_geo_info.csv")

# ----- drug data ------ #
# create the synonym DF
drug_full_info <- read.delim("data/db_data/drugbank_parsed.txt")
drug_info_df <- read.delim("data/db_data/drugbank_vocab_no_nutra.txt")
head(drug_info_df)

# only match drug data with more than 3 characters
drug_info_df$numchar <- sapply(drug_info_df$name, nchar)
short_names <- filter(drug_info_df, numchar <= 3)
names_to_match <- filter(drug_info_df, numchar > 3)


drug_info_df$nwords <- sapply(drug_info_df$name, function(x) length(strsplit(x, " ")[[1]]))

text_df <- gse_data[,c("gse", "str")]
text_df$str <- sapply(text_df$str, function(x)
{y <- gsub(";\t", " ", x);  ## fix ;\t
z <- gsub(' [[:punct:]]|[[:punct:]] ', ' ', y);  ##remove trailing punctuation
return(z)
})
gse_unigrams <- text_df %>% unnest_tokens(word, str, token=stringr::str_split, pattern = " ") 
gse_unigrams2 <- filter(gse_unigrams, str_length(word) > 3 & is.na(as.numeric(word))) # remove short words or numbers 
comb_names <- inner_join(filter(drug_info_df, nwords==1), gse_unigrams2, by=c("name"="word")) %>% distinct()
length(unique(comb_names$gse)) # 13608 GSEs
drugs <- unique(comb_names$name ) # 1309
gse_map_counts <- gse_unigrams %>% filter (word %in% drugs) %>% group_by(word) %>% summarise(total=n())  %>% arrange(desc(total))
View(gse_map_counts)

# TODO - decide if we are ok with discarding the list of stopwords from JAKE
#  problem - discarding aa from the drugbank classification
drug_stopwords <- c("same", "dmso", "water", "date", "sage", "atra", "camp", "saha", "balance", "date", "biotin")
jake_stopwords <- read.delim("data/ref_data/jake_stopwords.txt", head=FALSE)$V1
comb_names2 <- filter(select(comb_names, c("name", "dbID", "gse")), ! name  %in% union(drug_stopwords, jake_stopwords))
length(unique(comb_names2$gse)) # 9551
length(unique(comb_names2$name)) # 1250

# join with ATC
comb_names3 <- left_join(comb_names2, drug_full_info[,c("dbID", "ATC")], by="dbID")
write.csv(comb_names3, file="data/gse2drugbank.csv")


# MAY WANT TO FILTER: glucose, oxygen, nitrogen, calcium
