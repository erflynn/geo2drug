# code for constructing a matrix
# this is a wrapper for any matrix construction, and will be used by train_mat 
# and test_mat in the future

require('MetaIntegrator')
require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
idx <- as.numeric(args[2]) 
run_v <- args[3] # run identifier
dat.file <- args[4] # this is the file name with the list of GSEs
my.dir <- args[5] # parent directory to look in e.g. "data/03_silver_std/"

consensus_genes <- read.csv(sprintf("data/consensus_genes_%s.csv", organism))
consensus.genes <- sapply(consensus_genes$consensus.genes, as.character)

# // TODO - move this to a common location
extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}


# function to construct the matirx
loadGSEData <- function(gse, dat_info, my.dir, MIN.ROWS=10000){
  gse_dat <- dat_info[dat_info$gse==gse,]
  gse2 <- strsplit(gse, "-")[[1]][[1]]
  miceadds::load.Rdata(sprintf("%s/00_mat_files/%s_mat.RData", my.dir, gse2), "mat_obj")
  mat_obj2 <- mat_obj[[sprintf("%s_series_matrix.txt.gz", gse)]]
  gene_df <- apply(mat_obj2$gene_mat[,gse_dat$gsm], c(1,2), as.numeric)
  pheno_df <- mat_obj2$pheno[gse_dat$gsm,]
  class <- ifelse(gse_dat$text_sex=="female", 0, 1)
  names(class) <- gse_dat$gsm
  keys <- rownames(gene_df)
  names(keys) <- keys
  if (nrow(gene_df) < MIN.ROWS){
    print(sprintf("%s had too few rows, n=%s.", gse, nrow(gene_df)))
    return(NA)
  }
  
  # create a meta-object
  my.obj <- list("expr"=gene_df, "pheno"=pheno_df, "class"=class, "keys"=keys, "formattedName"=gse)
  if(!checkDataObject(my.obj, "Dataset")){
    print(sprintf("Error with %s", gse))
    return(NA)
  }
  
  # change the column labeling
  colnames(gene_df) <- 
    sapply(colnames(gene_df), 
           function(gsm) sprintf("%s.%s", gse, gsm))
  
  return(gene_df)
}


dat_info <- read_csv(dat.file)
list.gses <- dat_info$gse
chunk.gses <- extractChunk(list.gses)

# load all the data
loaded.gses <- lapply(chunk.gses,  function(gse) 
  loadGSEData(gse, dat_info, my.dir) )

# convert to expression matrix and write it out
gses_f <- loaded.gses[!is.na(loaded.gses)]
gses_expr <- lapply(gses_f, function(dat)
  exprsex::reorderRank(dat, gene_list=consensus.genes, to.rank=FALSE))
write_csv(gses_expr, sprintf("%s/03_out_mat/dat_%s_expr_%s.csv", my.dir, run_v, idx))

# extract and write out labels
dat_lab <- dat_info %>% filter(gse %in% chunk.gses) %>% mutate(sex=ifelse(text_sex=="male", 1, 0))
write_csv(dat_lab, sprintf("%s/03_out_mat/dat_%s_lab_%s.csv", my.dir, run_v, idx))

