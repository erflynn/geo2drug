# code for constructing the training and testing matrices

require('MetaIntegrator')
require('tidyverse')


##### ----- CONSTRUCT THE TRAIN MATRICES ----- ####

load("data/03_silver_std/03_out_mat/metaObj.RData") # --> metaObject

# iterate through the gses
train.gses <- unique(train_dat$gse)
trainIn <- metaObject$originalData
trainIn <- lapply(names(trainIn), function(gse) {
  x <- trainIn[[gse]];
  x$gene_mat <- x$expr;
  colnames(x$gene_mat) <- 
    sapply(colnames(x$gene_mat), 
           function(gsm) sprintf("%s.%s", str_replace_all(gse, "\\.", "-"), gsm))
  return(x)
  })
names(trainIn) <- names(metaObject$originalData)
consensus.genes <- exprsex::getConsensusGenes(list(trainIn)) # double nest for this
write_csv(data.frame(consensus.genes), "data/consensus_genes.csv")

train_expr <- lapply(trainIn, function(dat)
                      exprsex::reorderRank(dat$gene_mat, gene_list=consensus.genes, to.rank=FALSE))
train_rank <- lapply(trainIn, function(dat)
  exprsex::reorderRank(dat$gene_mat, gene_list=consensus.genes))

train_expr2 <- data.frame(do.call(cbind, train_expr))
train_rank2 <- data.frame(do.call(cbind, train_rank))

train_lab <- lapply(names(trainIn), function(gse) {
  x <- trainIn[[gse]];
  lab <- x$class;
  data.frame("gse"=rep(gse, length(lab)), "gsm"=names(lab), "sex"=lab)
  })
# output: gse, gsm, class label
train_lab2 <- do.call(rbind, train_lab)

write_csv(train_expr2, "data/03_silver_std/03_out_mat/train_expr.csv")
write_csv(train_rank2, "data/03_silver_std/03_out_mat/train_rank.csv")
write_csv(train_lab2, "data/03_silver_std/03_out_mat/train_lab.csv")


##### ----- CONSTRUCT THE TEST MATRICES ----- ####
test_dat <- read_csv("data/01_sample_lists/human_testing.csv")
MIN.ROWS <- 10000
test.gses <- unique(test_dat$gse)
test <- lapply(head(test.gses), 
                function(gse){
                  print(gse)
                  gse_dat <- test_dat[test_dat$gse==gse,]
                  gse2 <- strsplit(gse, "-")[[1]][[1]]
                  miceadds::load.Rdata(sprintf("data/03_silver_std/00_mat_files/%s_mat.RData", gse2), "mat_obj")
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
consensus_genes <- read_csv("data/consensus_genes.csv")
consensus.genes <- sapply(consensus_genes$consensus.genes, as.character)


test_expr <- lapply(test_f, function(dat)
  exprsex::reorderRank(dat, gene_list=consensus.genes, to.rank=FALSE))
test_rank <- lapply(test_f, function(dat)
  exprsex::reorderRank(dat, gene_list=consensus.genes))
test_expr2 <- data.frame(do.call(cbind, test_expr))
test_rank2 <- data.frame(do.call(cbind, test_rank))


write_csv(test_expr2, file="data/03_silver_std/03_out_mat/test_expr.csv")
write_csv(test_expr2, file="data/03_silver_std/03_out_mat/test_rank.csv")

filt_gses <- unique(sapply(colnames(text_expr2), function(x) strsplit(x, "\\.")[[1]][[1]]))
# extract test_lab from test_dat but filter out removed gses
test_lab <- test_dat %>% filter(gse %in% filt_gses) %>% mutate(sex=ifelse(sex=="male", 1, 0))
write_csv(test_lab, file="data/03_silver_std/03_out_mat/test_lab.csv")
