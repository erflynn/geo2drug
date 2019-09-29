# extract_gpl_consensus.R
# E Flynn
# 9/27/2019
#
# Code for extracting lists of GPLs and creating a consensus list of genes

source("code/utils/labeling_utils.R")



#' Extract a list of all GPLs and grab their data from GEOMetadb
#' @param df - the df to extract gpls from
#' @param df_id - the id to write results to
#' @con - connection to GEOMetadb
writeGPLs <- function(df, df_id, con){
  df_gpls <- dbGetQuery(con, sprintf("SELECT * FROM gpl WHERE gpl IN ('%s');",
                                     formattedList(df$gpl)))
  df_gpls %>% 
    select(gpl, title, manufacturer, technology, description, data_row_count) %>% 
    write_csv(sprintf("data/sample_lists/%s_gpl.csv", df_id))
  
}


# -- Make a consensus gene list using the genes in the top platforms -- #

#' Count the number of GSEs per GPL
#' 
#' @param df
#' @return gpls with the count of gses, arranged in descending order
countGPLs <- function(df) {
  df %>% 
    select(gse, gpl) %>% 
    unique() %>% 
    group_by(gpl) %>% 
    mutate(n=n()) %>% 
    select(-gse) %>% 
    unique() %>% 
    arrange(desc(n))
}



#' Extract genes from a list of platforms and filter by genes present
#' a pre-specified number of times to get a consensus gene list.
#' 
#' @param list.gpls - list of platforms to look at
#' @param cutoff.n - the minimum number of platforms a gene must be present in (default = 4)
#' 
#' @return the consensus gene list
consensusGenes <- function(list.gpls, cutoff.n=4){
  list.gpls <- hs_count_gpls$gpl[1:8]
  lists.genes <-  lapply(list.gpls, function(gpl) {
    keys <- exprsex:::.getGenes(gpl);
    gene.df <- exprsex:::.getGeneDf(keys); 
    filtered.df <- gene.df %>% select(gene) %>% unique();
    return(filtered.df)
  }) 
  df <- do.call(rbind,lists.genes) 
  print("Genes extracted, counting")
  gene_counts <- df %>% group_by(gene) %>% mutate(n=n()) %>% unique() %>% filter(n >= cutoff.n) 
  return(unique(gene_counts$gene))
}



# ------- MAIN ------ #

# read in data
mm2 <- read_csv("data/sample_lists/mouse_gse_gsm.csv")
rn2 <- read_csv("data/sample_lists/rat_gse_gsm.csv")
hs2 <- read_csv("data/sample_lists/human_gse_gsm.csv")

# extract and write out GPL info
con <- dbConnect(SQLite(), GEOMETADB_PATH)
writeGPLs(hs, "human", con)
writeGPLs(mm, "mouse", con)
writeGPLs(rn, "rat", con)
dbDisconnect(con)

# count gpls per gse
mm_count_gpls <- mm2 %>% countGPLs()
rn_count_gpls <- rn2 %>% countGPLs()
hs_count_gpls <- hs2 %>% countGPLs()

# get consensus genes for the top 8 platforms for each species
mm_gene_list <- consensusGenes(mm_count_gpls$gpl[1:8])
hs_gene_list <- consensusGenes(hs_count_gpls$gpl[1:8])
rn_gene_list <- consensusGenes(rn_count_gpls$gpl[1:8])

# write out consensus lists
write_csv(data.frame(mm_gene_list),  "data/ref_data/mm_consensus_list.csv", col_names=FALSE)
write_csv(data.frame(hs_gene_list),  "data/ref_data/hs_consensus_list.csv", col_names=FALSE)
write_csv(data.frame(rn_gene_list),  "data/ref_data/rn_consensus_list.csv", col_names=FALSE)

## TODO - there is something going wrong for rat and mouse :(


