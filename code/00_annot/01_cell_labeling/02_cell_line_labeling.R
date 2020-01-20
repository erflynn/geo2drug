# cell_line_labeling.R
# E Flynn
# 3/7/2019
# Updated: 4/12/2019
#
# Cell line labeling of ALL GSE data from human/mouse/rat.
#
# We are mapping if:
#  1) cell line is >= 3char and maps directly
#  2) if <= 3 char: 
#      THEN - must precede cell or cells or cell line
#  also considered cell line if description contains the phrase "cell line"
#
# Note that this mapping is at the GSE and not GSM level.

require('rjson')
require('tidytext')
require('tidyverse')
require('fuzzyjoin')
require('stringr')
require('Hmisc')
options(stringsAsFactors=FALSE)

# ---- GSE data ---- #
gse_data <- read.csv("data/01_sample_lists/gse_all_geo_info.csv")
gse_data$str <-sapply(1:nrow(gse_data), function(i) {paste(gse_data[i,c("title", "summary", "overall_design")][!is.na(gse_data[i,c("title", "summary", "overall_design")])], collapse=" ")}) 
gse_data$cell_line <- str_detect(gse_data$str, "cell line")
#gsm_filt_w_gse <- fread("data/00_db_data/gsm_all_geo_info.csv")

# ----- cell line data ------ #
# create the synonym DF
cell_info_df <- read.csv("data/00_db_data/cellosaurus_df.txt")
cell_lab <- cell_info_df %>% separate_rows(synonyms, sep="\\|") %>% select(synonyms, accession, cl)
cell_lab2 <- data.frame(apply(cell_lab, c(1,2) , function(x) str_trim(tolower(x))))
cell_lab_syn <- cell_lab2 %>% select(accession, synonyms) %>% filter(synonyms != "") 
cell_lab_name <- cell_lab2 %>% select(accession, cl) %>% distinct()
cell_df <- rbind(dplyr::rename(cell_lab_syn, cl=synonyms), cell_lab_name) 
cell_df$nwords <- sapply(cell_df$cl, function(x) length(strsplit(x, " ")[[1]]))

# divide into cell names >3 and less than 3
cell_df$numchar <- sapply(cell_df$cl, nchar)
cell_df2 <- cell_df[cell_df$numchar >3, ]
cell_df2_short <- cell_df[cell_df$numchar %in% c(2, 3), ]
cell_df2_short$cl <- sapply(cell_df2_short$cl, function(x) sprintf("%s cell", x))


# --- GSE data --- # 

# clean up text

text_df <- gse_data[,c("gse", "str")]
text_df$str <- sapply(text_df$str, function(x)
  {y <- gsub(";\t", " ", x);  ## fix ;\t
  z <- gsub(' [[:punct:]]|[[:punct:]] ', ' ', y);  ##remove trailing punctuation
  clean_str <- gsub('cells', 'cell', z); # convert cells --> cell
  return(clean_str)
})
gse_unigrams <- text_df %>% unnest_tokens(word, str, token=stringr::str_split, pattern = " ") 
gse_unigrams2 <- filter(gse_unigrams, str_length(word) > 3 & is.na(as.numeric(word))) # remove short words or numbers 
comb_names <- inner_join(filter(cell_df2, nwords==1), gse_unigrams2, by=c("cl"="word")) %>% distinct()

# filter for high TF
cells <- unique(comb_names$cl )
gse_map_counts <- gse_unigrams %>% filter (word %in% cells) %>% group_by(word) %>% summarise(total=n())  %>% arrange(desc(total))

# selected by looking at gse_map_counts
cl_stopwords <- c("cancer", "time", "center", "rare", "focus", "peak", "sage", "bona", "mast", "fisher", "bones", "patches", "madison", "ears", "chance", "cost")
comb_names2 <- comb_names %>% filter(!cl %in% cl_stopwords)
length(unique(comb_names$gse)) # 11251 out of 45036

comb_data <- right_join(comb_names2[,c("accession", "cl", "gse")], gse_data[,c("gse", "cell_line")])
write.table(comb_names2[,c("accession", "cl", "gse")], file="data/tmp_one_word_cell_mapping.txt", row.names=FALSE, sep="\t")


# NOW SEARCH FOR LONGER STRINGS 
cell_long <- cell_df2[cell_df2$nwords > 1,]
cell_long2 <- rbind(cell_long, cell_df2_short)

text_df$str <- gsub("\t", " ", text_df$str)
text_df$str <- gsub('"', "", text_df$str)

text_df2 <- text_df
text_df2$str <- escapeRegex(text_df$str)
cell_long2$cl <- sapply(escapeRegex(cell_long2$cl), function(x) sprintf(" %s ", str_trim(x)))

# write these out - will join by exact match in python (so slow :,( )
write.table(cell_long2[,c("accession", "cl")], file="data/tmp/multiword_cell.txt", row.names=FALSE, sep="\t")
write.table(text_df2, file="data/tmp_gse_text.txt", row.names=FALSE, sep="\t")

# bigrams, trigrams
bigrams <- gse_data[,c("gse", "str")] %>% unnest_tokens(bigram, str, token = "ngrams",  n = 2)
bigrams_separated <- bigrams %>% separate(bigram, c("word1", "word2"), sep = " ")
expanded_stop_words <- c(stop_words$word, "human", "hospital")
bigrams_filt <- bigrams_separated %>% 
  filter(is.na(as.numeric(word1)))  %>% filter(is.na(as.numeric(word2))) %>%
  filter(!word1 %in% expanded_stop_words) %>%
  filter(!word2 %in% expanded_stop_words)
# filter for stopwords, numeric, etc
bigram_counts <- bigrams_filt %>% 
  count(word1, word2, sort = TRUE)

bigrams_united <- bigrams_filt %>%
  unite(bigram, word1, word2, sep = " ")


cell_two_words <- rbind(filter(cell_df2, nwords==2), cell_df2_short)
cells2 <- unique(cell_two_words$cl)
bigram_map_counts <- bigrams_united %>% filter (bigram %in% cells2) %>% group_by(bigram) %>% summarise(total=n())  %>% arrange(desc(total))
mapped_two <- inner_join(bigrams_united, cell_two_words, c("bigram"="cl"))
# "apoe ko", "renal carcinoma", "b16 melanoma"


## three ##
trigrams <- gse_data[,c("gse", "str")] %>% unnest_tokens(trigram, str, token = "ngrams", n = 3)
cell_tri <- filter(cell_df2, nwords>=3) %>% unnest_tokens(cl, cl, token="ngrams", n=3)
cell_three_words <- rbind(filter(cell_df2, nwords==3), cell_df2_short)
mapped_three <- inner_join(trigrams, cell_three_words, c("trigram"="cl"))

trigram_counts <- mapped_three %>% group_by(trigram) %>% summarise(total=n())  %>% arrange(desc(total))


### put together all the mapping
mapped_dat <- do.call(rbind, list(comb_names2[,c("accession", "cl", "gse")],
  mapped_two %>% select(gse, bigram, accession) %>% dplyr::rename(cl=bigram),
  mapped_three %>% select(gse, trigram, accession) %>% dplyr::rename(cl=trigram)))


mapped_all <- right_join(mapped_dat, gse_data[,c("gse", "cell_line")])

df <- data.frame(cbind(mapped_all$cell_line, !is.na(mapped_all$cl)))
colnames(df) <- c("cell line", "mapped")
table(df)

# write out
write.table(mapped_all, file="data/02_labeled_data/cell_line_mapped_gse.txt", row.names=FALSE)
gse_to_keep <- filter(mapped_all, !cell_line & is.na(cl))$gse
write.table(data.frame(gse_to_keep), file="data/02_labeled_data/non_cell_line_gse.txt", row.names=FALSE)
