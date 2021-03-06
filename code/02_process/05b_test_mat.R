# code for constructing the training and testing matrices

require('MetaIntegrator')
require('tidyverse')

SIZE.CHUNK <- 30
args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
idx <- as.numeric(args[2])
run_v <- args[3]

consensus_genes <- read.csv(sprintf("data/consensus_genes_%s.csv", organism))
consensus.genes <- sapply(consensus_genes$consensus.genes, as.character)

if (run_v == "common"){
  test_dat <- read_csv(sprintf("data/01_sample_lists/%s_testing.csv", organism))
} else {
  test_dat <- read_csv(sprintf("data/01_sample_lists/%s_testing_%s.csv", organism, run_v))
}


test.gses <- unique(test_dat$gse)
NUM.CHUNKS <- ceiling(length(test.gses)/SIZE.CHUNK)
end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(test.gses), (idx+1)*SIZE.CHUNK)

test.gses2 <- test.gses[(idx*SIZE.CHUNK+1):end_idx]



##### ----- CONSTRUCT THE TEST MATRICES ----- ####

MIN.ROWS <- 10000

# divide into groups of 20
test <- lapply(test.gses2, 
                function(gse){
                  print(gse)
                  gse_dat <- test_dat[test_dat$gse==gse,]
                  gse2 <- strsplit(gse, "-")[[1]][[1]]
                  miceadds::load.Rdata(sprintf("data/03_silver_std/%s/00_mat_files/%s_mat.RData", organism, gse2), "mat_obj")
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
                })

test_f <-test[!is.na(test)]


test_expr <- lapply(test_f, function(dat)
  exprsex::reorderRank(dat, gene_list=consensus.genes, to.rank=FALSE))
test_rank <- lapply(test_f, function(dat)
  exprsex::reorderRank(dat, gene_list=consensus.genes))
save(test_expr, file=sprintf("data/03_silver_std/%s/03_out_mat/test_%s_expr_%s.RData", organism, run_v, idx))
save(test_rank, file=sprintf("data/03_silver_std/%s/03_out_mat/test_%s_rank_%s.RData", organism, run_v, idx))
test_expr2 <- data.frame(do.call(cbind, test_expr))
test_rank2 <- data.frame(do.call(cbind, test_rank))


write_csv(test_expr2, sprintf("data/03_silver_std/%s/03_out_mat/test_%s_expr_%s.csv", organism, run_v, idx))
write_csv(test_rank2, sprintf("data/03_silver_std/%s/03_out_mat/test_%s_rank_%s.csv", organism, run_v, idx))

filt_gses <- unique(sapply(colnames(test_expr2), function(x) strsplit(x, "\\.")[[1]][[1]]))
# extract test_lab from test_dat but filter out removed gses
test_lab <- test_dat %>% filter(gse %in% filt_gses) %>% mutate(sex=ifelse(text_sex=="male", 1, 0))
write_csv(test_lab, sprintf("data/03_silver_std/%s/03_out_mat/test_%s_lab_%s.csv", organism, run_v, idx))
