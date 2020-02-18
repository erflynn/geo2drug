# Code for extracting consensus genes

require('MetaIntegrator')
require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]

# consensus genes for seq data
load(sprintf("../rnaseq_data/%s_matrix.rda", organism))
gene_meta <- meta[,c("genes", "gene_entrezid")] # the genes are the rownames for the expression mat
# // TODO: do we want the "reads_aligned" field? do we want to normalize here?



# consensus genes for platform data
# for each of the top ten platforms, read in one study per platform
grp.df <- read_csv(sprintf("data/01_sample_lists/%s_oligo_grps.csv", organism))
grps <- grp.df %>% filter(grp %in% c("rep_sample", "mixed_sex")) %>% group_by(gpl) %>% 
  select(gse) %>% group_split()

getGSE <- function(gse){
  my.f <- sprintf("data/10_oligo/%s/00_mat_files/%s_mat.RData", organism, gse)
  if (!file.exists(my.f)){
    return(NA)
  }
  miceadds::load.Rdata(
    my.f,
    "gse_mat") 
  if (length(gse_mat) != 1){
    return(NA)
  }
  if (is.na(gse_mat[[1]])){
    return(NA)
  }
  genes <- rownames(gse_mat[[1]]$gene_mat)
  print(length(genes))
  if (length(genes) < 8000){
    return(NA)
  } 
  return(gse_mat[[1]])
}


gse.mats <- sapply(grps,
function(grp) {
  gpl.gses <- grp$gse
  print(gpl.gses)
  my_mat <- NA
  i <- 1
  while(is.na(my_mat) & i <= length(gpl.gses)){
    gse <- gpl.gses[[i]]
    print(gse)
    my_mat <- getGSE(gse)
    i = i+1
  }
  return(my_mat)
})


#### DEPRECATED ####

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
