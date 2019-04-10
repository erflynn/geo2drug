# cell_line_labeling.R
# E Flynn
# 3/7/2019
#
# Cell line labeling

require('rjson')

# load the cellosaurus labels
cellosaurus <- fromJSON(file="data/db_data/cellosaurus.json") 
cell_names <- names(cellosaurus)
names2 <- unlist(sapply(cellosaurus, function(x) return(x$name)))
cell_syns <- unlist(sapply(cellosaurus, function(x) return(x$synonyms)))
cell_lst <- c(cell_names, cell_syns)

# data frame
cl_df <- data.frame("cl"=cell_names, "syn"=sapply(cellosaurus, function(x) paste(x$synonyms, collapse=" ; ")))
cl_df <- separate_rows(cl_df, syn, sep=";")

# now string trim everything
cl_df$syn <- sapply(cl_df$syn, str_trim)
cl_df$cl <- sapply(cl_df$cl, str_trim)



# load the ALE labels - does this help?
ale_data <- read.csv("../../drug_expression/drug_labeling/ale_processing/ale_combined_data.csv", stringsAsFactors = FALSE)
cell_line_intersect <- intersect(sapply(ale_data$text_tissue_name, tolower), sapply(cell_lst, tolower)) # ok... not useful

# map text tissue to cellosaurus ID
ale_tiss <- data.frame("name" =unique(ale_data[,c("text_tissue_name")]), "tokens"=unique(ale_data[,c("text_tissue_name")]))
ale_tiss <- separate_rows(ale_tiss, tokens, sep=" ")
ale_to_cl <- inner_join(ale_tiss, cl_df, by=c("tokens"="cl"))
ale_to_cl2 <- inner_join(ale_tiss, cl_df, by=c("tokens"="syn"))
ale_to_cl <- ale_to_cl[, c("name", "tokens")]
mapped1 <- ale_to_cl[!duplicated(ale_to_cl),]
mapped1 <- rename(mapped1, "cl"="tokens")
ale_to_cl2 <- ale_to_cl2[, c("name", "cl")]
mapped2 <- ale_to_cl2[!duplicated(ale_to_cl2),]

text_to_cl <- rbind(mapped1, mapped2)
text_to_cl <- text_to_cl[!duplicated(text_to_cl),]
text_to_cl <- text_to_cl[text_to_cl$name != "",] # 311, this removes a ton of junk
length(unique(ale_data$text_tissue_name)) # 1594
length(unique(text_to_cl$name)) #309

ale_cell_line <- ale_data[,c("text_tissue_name", "cell_line")]
ale_cell_line <- ale_cell_line[!duplicated(ale_cell_line),]
text_to_cl <- rename(text_to_cl, "text_tissue_name"="name")
ale_cl_mapping <- left_join(ale_cell_line, text_to_cl)
table(ale_cl_mapping$cell_line, is.na(ale_cl_mapping$cl))
unmapped_cell_line <- filter(ale_cl_mapping, is.na(cl) & cell_line==TRUE) # cell line that don't map to cellosaurus 
new_cell_line <- filter(ale_cl_mapping, !is.na(cl) & cell_line==FALSE) # cell line that do map but weren't listed as cell line
not_cell_line <- filter(ale_cl_mapping, is.na(cl) & cell_line==FALSE) # 963
cell_word <- sapply(not_cell_line$text_tissue_name, function(x) "cell" %in% strsplit(x, " ")[[1]])
not_cell_line_hc <- not_cell_line[!cell_word,] 

possible_map <- not_cell_line[cell_word,] # M-07E, Th1, Th2, HeLa-S3, M, MM1-S, 22Rv-1, A-375P, Th17, etc...
#  TODO
#  - some of these have a substring that *is* a cell line name
# filter(ale_cl_mapping, !is.na(cl) & cell_line==TRUE) # <-- HC mapping
not_cell_line
head(cell_line)
head(text_to_cl)



# what fraction of these are actually in the same data 

cell_info_df <- do.call(rbind, 
                   lapply(cellosaurus, function(x) 
                     data.frame(lapply(x, function(y) paste(y, collapse=" | ")), stringsAsFactors=FALSE)))
# add the sex labels in 
cell_info_df$accession <- sapply(cell_info_df$accession, function(x) strsplit(x, " \\| ")[[1]][[1]])

cell_info_df$cl <- names(cellosaurus)
write.csv(cell_info_df, "data/db_data/cellosaurus_df.txt", row.names=FALSE)

tiss_names <- sapply(unique(ale_data$text_tissue_name), function(x) strsplit(x, " ", fixed=TRUE)[[1]])
overlap <- intersect(unlist(tiss_names), cell_lst) # 450 overlap!!!

# map the ALE data to a cellosaurus ID
# gse ale_text_tissue cellosaurus_name cellosaurus_id
mapTextTissToCL <- function(text_name){
  name_tokens <- strsplit(text_name, " ", fixed=TRUE)[[1]]
  cl_names <- intersect(name_tokens, overlap)
  if (length(cl_names)==0){
    return(NA)
  }
  cl_name <- intersect(cl_df$cl, cl_names)
  if (length(cl_name)!=0){
    return(cl_name)
  }
  syn_name <- filter(cl_df, syn %in% cl_names)$cl
  if (length(syn_name) != 0){
    return(syn_name)
  }
  return(NA)
}

cl_labels <- sapply(unique(ale_data$text_tissue_name), mapTextTissToCL) # this is slow... speed up
table(is.na(cl_labels))
names(cl_labels) <- ale_data$text_tissue_name
cl_labels2 <- cl_labels[!is.na(cl_labels)]
cl_labels3 <- (cl_labels2[!duplicated(cl_labels2)]) # 447
table(sapply(cl_labels3, length))
# disambiguate these two, discard the rest
cl_labels3$`U-373MG cell` <- "U-373MG ATCC"
cl_labels3$`U-87MG cell` <- "U-87MG ATCC"

cl_labels4 <- cl_labels3[sapply(cl_labels3, length) == 1]
cl_df_ale <- data.frame(cbind("cell"=cl_labels4, "text_tissue_name"=names(cl_labels4)))
head(ale_cl_mapping)


# COMBINE THIS WITH THE OTHER DATA




#unique(filter(ale_data, cell_line==TRUE)$text_tissue_name)

# what is the count by sex
comb <- inner_join(ale_cl_mapping, cell_info_df)
dim(comb)
table(comb$sex)
write.table(comb, file="data/rough_cl_mapping.txt", sep="\t", row.names=FALSE)

# load the MeSH IDs


require('fuzzyjoin')
# not sure I actually want to use this... 
# -- problems - inexact matches with cell lines
# -- I don't think it's worth it to be doing this
# also a problem:
#   we're mapping: BTO --> text --> cell line
#    this means that if the cell line isnt in BTO, then it doesnt work 


# try mapping the manual data - what cell lines do we see?
# the problem is that these were also labeled w BTO
manual_labels <- read.csv("../tissue_labeling/data/manual_w_category.csv")
map_manual<- sapply(manual_labels$tissue_name, mapTextTissToCL)
manual_labels$cell_line <- map_manual
table(is.na(map_manual))
table(manual_labels[!is.na(map_manual),]$tissCat) # all of these were mapped to "Other"
manual_labels$contains_cell <- sapply(manual_labels$tissue_name, function(x) "cell" %in% strsplit(x, " ")[[1]])
table(!is.na(manual_labels$cell_line), manual_labels$"contains_cell") 
head(unique(manual_labels$cell_line),50)
# 2787 cell line
# no examples fall into cell line and doesn't contain the word cell
unique(manual_labels[is.na(manual_labels$cell_line) & manual_labels$"contains_cell",]$tissue_name) # probs mislabeled?
#BALL-1, TALL-1, IOSE, HBMEC, KOPT-K1 <-- these should map I think, but they dont

num_terms <- sapply(cl_df$cl, function(x) length(strsplit(x, " ")[[1]]))

# read in the tissue states manual labels
tiss_labels <- read.csv("../tissue_labeling/data/manual_acc_0119.csv")

gses.remove <- read.csv("../tissue_labeling/data/gses_to_remove.csv")
tiss_labels_filt <- tiss_labels[!tiss_labels$gse %in% gses.remove[,1],]
head(tiss_labels_filt)
tiss_names <- c("Adipose", "Blood", "Brain", "Breast", "Colon", 
"Heart", "Kidney", "Liver", "Lung","Muscle", "Prostate", "Skin")
getMaxTissue <- function(ts, tiss_names){
  apply(ts[,tiss_names], 1, function(y) tiss_names[unlist(which.max(y))])
}
tiss_labels_filt$maxTiss <- getMaxTissue(tiss_labels_filt, tiss_names)
manual_labels <- rename(manual_labels, "gsm"="X")
comb_manual <- full_join(manual_labels, select(tiss_label, "gsm", "gse", "maxTiss"))
tmp <- filter(comb_manual, ((!is.na(cell_line)) ))
tmp2 <- tmp[sapply(tmp$maxTiss, function(x) !is.null(x)),] # why did so many of these not run? :/ - they also have no GSE listed
table(tmp$cell_line) # there are 35 cell lines - we only have data from 5 of them - what happened here?
tmp3 <- tmp2[sapply(tmp2$maxTiss, function(x) length(x)!=0),]
head(tmp3)
table(tmp3$cell_line, sapply(tmp3$maxTiss, unlist))
#           Blood Brain Breast Kidney Skin
# CCRF-CEM     5     0      0      0    0 # T-ALL
# MOLT-4       5     0      0      0    0 # T-ALL
# SaOS-2       0     0      6      0    0 # osteosarcoma
# SK-OV-3      0     0     10      1    0 # ovarian
# U2OS         5    22     24      0    1 # osteosarcoma
counts.per.cat <- table(apply(tmp2[,c("tissue_name", "maxTiss")], 1, function(x) paste(x, collapse="\t")))
counts.per.cat[order(counts.per.cat)]


# SEX LABELS FOR CELL LINES - what is the breakdown/distribution?

