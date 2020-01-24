# the goal of this is to assess how much missingness is ok

# look at the the training data from the test platforms



require('exprsex')
require('tidyverse')
require('data.table')
source("code/utils/assessment_utils.R")

# <--- READ IN THE INPUT ---> #
options(stringsAsFactors = FALSE)

# ---- prepare data --- #
common.genes <- read.csv(sprintf("data/consensus_genes_%s.csv", organism))$consensus.genes
if (organism=="human"){
  full.genes <- read.csv(sprintf("data/consensus_genes_%s_full.csv", organism))$consensus.genes
  genes <- full.genes
} else {
  genes <- common.genes
}


train_rank <- read_process_mat(sprintf("data/03_silver_std/%s/03_out_mat/train_%s_rank.csv", organism, run_v), genes)
train_expr <- read_process_mat(sprintf("data/03_silver_std/%s/03_out_mat/train_%s_expr.csv", organism, run_v), genes)
train_lab <- read_process_lab(sprintf("data/03_silver_std/%s/03_out_mat/train_%s_lab.csv", organism, run_v))



sl_genes <- read_clean_sl_genes(organism, run_v)

# ok now we want the counts of missingness by column
sl_overlap <-intersect(c(sl_genes$f, sl_genes$m), rownames(train_rank))
train_sl <- train_rank[sl_overlap,] # there are 31 total
nas.by.col <- apply(train_sl, 2, function(x) sum(is.na(x)))
summary(nas.by.col)
table(nas.by.col<=1) # 4510 <-- we can try these!
load(file=sprintf("data/04_fits/fit_%s_%s.RData", organism, run_v)) # full
preds.1na <- predSexLab(fit, as.matrix(train_sl[,c(nas.by.col<=1)]))
calcAcc(train_lab[c(nas.by.col<=1)], preds.1na) # 96.9%

# missingness by platform??

dat <- train_sl[,c(nas.by.col<=1)]
test_labels <- train_lab[c(nas.by.col<=1)]

set.seed(116)
removeGenes <-function(nGenes, dat, test_labels) {
  to.remove <- sl_overlap[sample(1:length(sl_overlap), nGenes, replace=FALSE)]
  dat[to.remove,] <- NA
  acc <- calcAcc(test_labels, predSexLab(fit, as.matrix(dat)))
  return(data.frame("ngenes"=nGenes, "acc"=acc))
}

## // TODO -- do this with the testing data
acc_remove <- do.call(rbind, lapply(1:50, function(i)
                     do.call(rbind, 
                      lapply(1:length(sl_overlap), function(x) removeGenes(x, dat, test_labels)))))

# change NAs to zero
acc_remove[is.na(acc_remove)] <- 0

acc_remove %>% write_csv("data/05_performance/human_missing_train.csv")
# PLOT
ggplot(acc_remove, aes(y=acc, x=factor(ngenes)))+geom_boxplot()+ylab("accuracy (50 runs)")+
  ylim(0,1.0)+xlab("number of genes removed")+
  ggtitle("Sensitivity to gene dropout")
ggsave(sprintf("figures/dropout_sens_%s_train.png", organism), dpi="print", width=6, height=4)

# look at the test data
list.test_c <- lapply(0:2, function(i) readInTest(organism, i, "common", common.genes)) 
test_rank <- do.call(cbind, lapply(list.test_c, function(x) x$rank))
test_lab <- unlist(sapply(list.test_c, function(x) x$lab))

test_sl <- test_rank[sl_overlap,] # there are 31 total
nas.by.col <- apply(test_sl, 2, function(x) sum(is.na(x)))
summary(nas.by.col)
table(nas.by.col<=1) # 4510 <-- we can try these!

test_rank2 <- test_sl[,c(nas.by.col<=1)]
test_lab2 <-test_lab[c(nas.by.col<=1)]

preds.1na.test <- predSexLab(fit, as.matrix(test_rank2))
calcAcc(test_lab2, preds.1na.test) # 95.3%

acc_remove_test <- do.call(rbind, lapply(1:50, function(i)
  do.call(rbind, 
          lapply(1:length(sl_overlap), function(x) removeGenes(x, test_rank2, test_lab2)))))

# change NAs to zero
acc_remove_test[is.na(acc_remove_test)] <- 0
acc_remove_test %>% write_csv("data/05_performance/human_missing_test.csv")

# PLOT
ggplot(acc_remove_test, aes(y=acc, x=factor(ngenes)))+geom_boxplot()+ylab("accuracy (50 runs)")+
  ylim(0,1.0)+xlab("number of genes removed")+
  ggtitle("Sensitivity to gene dropout")
ggsave(sprintf("figures/dropout_sens_%s_test.png", organism), dpi="print", width=6, height=4)


# // TODO: 
# - accuracy vs number of genes missing for the test data?
# - gene importance

