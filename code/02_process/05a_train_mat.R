# code for constructing the training and testing matrices

require('MetaIntegrator')
require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
run_v <- args[2]

##### ----- CONSTRUCT THE TRAIN MATRICES ----- ####
train_dat <- read_csv(sprintf("data/01_sample_lists/%s_training.csv", organism))
if (run_v == "full"){
  train_dat2 <- read_csv(sprintf("data/01_sample_lists/%s_training_full.csv", organism))
  train_dat <- rbind(train_dat, train_dat2)
}
miceadds::load.Rdata(sprintf("data/03_silver_std/%s/04_meta_res/input_metaObj_%s.RData", organism, run_v), "metaObj") 

# iterate through the gses
train.gses <- unique(train_dat$gse)
trainIn <- metaObj$originalData
trainIn <- lapply(names(trainIn), function(gse) {
  x <- trainIn[[gse]];
  x$gene_mat <- x$expr;
  colnames(x$gene_mat) <- 
    sapply(colnames(x$gene_mat), 
           function(gsm) sprintf("%s.%s", str_replace_all(gse, "\\.", "-"), gsm))
  return(x)
  })
names(trainIn) <- names(metaObj$originalData)
# I want to use the same set of consensus genes for EVERYTHING
# this means I will define the consensus genes based on the "full" list and just leave it as such 
# for human data, for the others -- I'll just use whatever they get
if (run_v == "full" | organism != "human"){
  consensus.genes <- exprsex::getConsensusGenes(list(trainIn)) # double nest for this
  write_csv(data.frame(consensus.genes), sprintf("data/consensus_genes_%s.csv", organism))
} else {
  consensus_genes <- read_csv(sprintf("data/consensus_genes_%s.csv", organism))
  consensus.genes <- sapply(consensus_genes$consensus.genes, as.character)
}

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

# we have to use write.csv to keep the rownames -- can grab them later if we want from
# consensus genes but it is something to be aware of
write.csv(train_expr2, file=sprintf("data/03_silver_std/%s/03_out_mat/train_%s_expr.csv", organism, run_v))
write.csv(train_rank2, file=sprintf("data/03_silver_std/%s/03_out_mat/train_%s_rank.csv", organism, run_v))
write_csv(train_lab2, sprintf("data/03_silver_std/%s/03_out_mat/train_%s_lab.csv", organism, run_v))

