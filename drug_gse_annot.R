# drug_gse_annot.R
# E Flynn
# 3/3/2019
#
# look at the breakdown of the annotations
#  - what is the sex breakdown? how many have sex labels?
#     - of these, what is the pattern?
#  - what is the tissue breakdown

require('tidyverse')
require('GEOmetadb')
options(stringsAsFactors=FALSE)


# --- load relevant data --- #

# text-based labels
ale_data <- read.csv("../../drug_expression/drug_labeling/ale_processing/ale_combined_data.csv")

# gse to mesh data
mapping_tab6 <-read.table("data/mesh_db_mapping_0302.txt", header=TRUE)
gse_mesh <- read.delim("data/gse_to_mesh.txt", header=TRUE)
gse_mesh_db <- full_join(gse_mesh, mapping_tab6, by="MeSH")
head(gse_mesh_db)

# expression-based labels
label_mat <- read.csv("data/label_mat.csv")
label_mat <- rename(label_mat, "gsm" ="X")

# combine
comb_labels <- left_join(label_mat, select(ale_data, c("gsm", "gse", "gpl", "text_sex", "text_tissue_name", "cell_line")))          
length(unique(comb_labels$gse)) # 1993
table(comb_labels$rank_sex) # 100,831
# female   male 
# 49795  51036 
table(comb_labels$toker_sex) # 59,531 (41,300 could not label) - 1,177 studies
# female   male 
# 28526  31005 
length(unique(filter(comb_labels, !is.na(toker_sex))$gse))

table(comb_labels$massir_sex) # 98,576 (2,255 could not label) - 1,907 studies
# female   male 
# 50707  47869 

length(unique(filter(comb_labels, !is.na(massir_sex))$gse)) # 35,659 samples (44,149) - 613 studies
table(comb_labels$text_sex) 
#         F     M 
# 44149 17062 18597 

# --- Make a bar chart of the sex breakdown --- #
sex_lab_present <- filter(comb_labels, text_sex %in% c("F", "M"))
sex_lab_present$text_sex <- ifelse(sex_lab_present$text_sex=="M", "male", "female") # fix this
length(unique(sex_lab_present$gse)) # 613 (1380 did not label)

# male only, female only, mixed 

num_studies <- length(unique(comb_labels$gse))

getSexLabStats <- function(df){
  # function to get stats of numbers m,f for a particular study 
  count.per.study <- df %>% group_by(gse) %>% summarize(num_f=sum(sexlab=="female"), num_m=sum(sexlab=="male"))
  m.only <- filter(count.per.study, num_f==0)
  f.only <- filter(count.per.study, num_m==0)
  mixed.sex <- filter(count.per.study, num_m!=0, num_f!=0) 
  
  num_studies <- length(unique(comb_labels$gse))
  list.stats <- list("f"=nrow(f.only), "m"=nrow(m.only), "both"=nrow(mixed.sex), "unlab"=num_studies-length(unique(df$gse)))
  return(list.stats)
}
sex_lab_present2 <- sex_lab_present
sex_lab_present2 <- rename(sex_lab_present2, sexlab=text_sex)
text.stats <- getSexLabStats(sex_lab_present2)
toker <- filter(comb_labels, !is.na(toker_sex))
toker <- rename(toker, sexlab=toker_sex)
toker.stats <- getSexLabStats(toker)
massir <- filter(comb_labels, !is.na(massir_sex))
massir <- rename(massir, sexlab=massir_sex)
massir.stats <- getSexLabStats(massir)
rank_dat <- rename(comb_labels, sexlab=rank_sex)
rank.stats <- getSexLabStats(rank_dat)

stat.table <- do.call(rbind, lapply(list(text.stats, toker.stats, massir.stats, rank.stats), data.frame))
stat.table$method <- c("text", "toker", "massir", "rankS")

melted.tab <- melt(stat.table)
head(melted.tab)
melted.tab$method <- factor(melted.tab$method, c("text", "toker", "massir", "rankS"))
ggplot(melted.tab, aes(x=method, y=value, fill=variable))+geom_bar(stat="identity")+ylab("Number of studies")+theme(text = element_text(size=20))



# ---- Calculate the accuracy ---- #

table(sex_lab_present$text_sex==sex_lab_present$toker_sex)  
# FALSE  TRUE 
# 3713 17438 --> 48.9% accuracy (82.4% of labeled accurate)
table(sex_lab_present$text_sex==sex_lab_present$massir_sex)
# FALSE  TRUE 
# 8370 26750 --> 75.0% accuracy (76.1% of labeled accurate)
table(sex_lab_present$text_sex==sex_lab_present$rank_sex) 
# FALSE  TRUE 
# 5091 30054  --> 84.2% accuracy

# look at concordance and mismatched data
table(comb_labels$toker_sex==comb_labels$massir_sex) # 47069 (80% concordance)
table(comb_labels$toker_sex==comb_labels$rank_sex) # 45711 (76.8% concordance)
table(comb_labels$massir_sex==comb_labels$rank_sex) # 71960 (73.2% concordance)

mismatched <- filter(sex_lab_present, text_sex!=rank_sex)
table(mismatched$massir_sex==mismatched$rank_sex) # 3069/5206 = 61.1%
table(mismatched$toker_sex==mismatched$rank_sex) # 1214/2145 = 56.5%


# ---- Plot tissue break down of data ---- #
tiss_dat <- comb_labels %>% group_by(gse) %>% summarize(common_tiss=names(which.max(table(tissue))))
ggplot(tiss_dat, aes(x=common_tiss, fill=common_tiss))+geom_histogram(stat="count")+xlab("Tissue breakdown")+ylab("Number of studies")+theme(text = element_text(size=20))

# what about for single sex?
count.per.study <- comb_labels %>% group_by(gse) %>% summarize(num_f=sum(rank_sex=="female"), num_m=sum(rank_sex=="male"))
m.only <- filter(count.per.study, num_f==0)
f.only <- filter(count.per.study, num_m==0)
mixed.sex <- filter(count.per.study, num_m!=0, num_f!=0) 

count.per.study$study_type <- ifelse(count.per.study$num_f==0, "m", ifelse(count.per.study$num_m==0, "f", "both"))
study_tiss_count <- left_join(count.per.study, tiss_dat)

# quick cell line filter
study_cell_line <- comb_labels %>% group_by(gse) %>% summarize(cell_line=any(cell_line))
table(study_cell_line$cell_line)
non_cell_line <- filter(study_cell_line, cell_line != "TRUE")$gse # this is not v accurate --> follow up w labels

# group by study type 
study_tab <-study_tiss_count[,c("gse", "study_type", "common_tiss")]
count_table <- table(study_tab[,c("study_type", "common_tiss")])
chisq.p <- lapply(colnames(count_table), function(x)
       chisq.test(count_table[c("f", "m"),x])$p.value)
names(chisq.p) <- colnames(count_table)
chisq.p[chisq.p < 0.05/ncol(count_table)]

ggplot(filter(study_tiss_count, !is.na(study_type) & study_tiss_count$gse %in% non_cell_line), aes(x=common_tiss, fill=common_tiss))+
  geom_histogram(stat="count")+xlab("Tissue breakdown")+ylab("Number of studies")+theme(text = element_text(size=20))+facet_grid(study_type ~ .)

# look at a couple problematic studies
left_join(filter(study_tiss_count, common_tiss=="Prostate" & study_type=="f" & gse %in% non_cell_line), comb_labels)
left_join(filter(study_tiss_count, common_tiss=="Breast" & study_type=="m"), comb_labels)


