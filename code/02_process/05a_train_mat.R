# code for constructing the training and testing matrices

require('MetaIntegrator')
require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]

##### ----- CONSTRUCT THE TRAIN MATRICES ----- ####
train_dat <- read_csv(sprintf("data/01_sample_lists/%s_training.csv", organism))
load(sprintf("data/03_silver_std/%s/03_out_mat/metaObj.RData", organism)) # --> metaObject

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
write_csv(data.frame(consensus.genes), sprintf("data/consensus_genes_%s.csv", organism))
#consensus_genes <- read.csv("data/consensus_genes.csv")
#consensus.genes <- sapply(consensus_genes$consensus.genes, as.character)


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

write_csv(train_expr2, sprintf("data/03_silver_std/%s/03_out_mat/train_expr.csv", organism))
write_csv(train_rank2, sprintf("data/03_silver_std/%s/03_out_mat/train_rank.csv", organism))
write_csv(train_lab2, sprintf("data/03_silver_std/%s/03_out_mat/train_lab.csv", organism))

