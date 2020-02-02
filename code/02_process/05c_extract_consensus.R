# Code for extracting consensus genes

require('MetaIntegrator')
require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]

##### ----- CONSTRUCT THE TRAIN MATRICES ----- ####
train_dat <- read_csv(sprintf("data/01_sample_lists/%s_training.csv", organism))
if (organism == "human"){ # for human, we use the full data, for the rest we use the common/smaller dataset
  train_dat2 <- read_csv(sprintf("data/01_sample_lists/%s_training_full.csv", organism))
  train_dat <- rbind(train_dat, train_dat2)
  miceadds::load.Rdata(sprintf("data/03_silver_std/%s/04_meta_res/input_metaObj_full.RData",organism), "metaObj") 
} else {
  miceadds::load.Rdata(sprintf("data/03_silver_std/%s/04_meta_res/input_metaObj_common.RData",organism), "metaObj") 
}

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

consensus.genes <- exprsex::getConsensusGenes(list(trainIn)) # double nest for this
write_csv(data.frame(consensus.genes), sprintf("data/consensus_genes_%s.csv", organism))
