# the goal of this is to assess how much missingness is ok

# look at the the training data from the test platforms

require('exprsex')
require('tidyverse')
require('data.table')
source("code/utils/assessment_utils.R")

# <--- READ IN THE INPUT ---> #
options(stringsAsFactors = FALSE)

# ---- prepare data --- #

ss <- read_process_ds(organism, "06_single_sex", "single_sex")
train <- read_process_ds(organism, "03_silver_std", "train_common")
train_full <- read_process_ds(organism, "03_silver_std", "train_full")
test <- read_process_ds(organism, "03_silver_std", "testing")
test_full <- read_process_ds(organism, "03_silver_std", "testing_full")


# ---- pre-process ---- #
train_rank <- exprsex:::.expDataToRanks(train[[1]]$expr)
train_lab <- train[[1]]$lab

test1_rank <- exprsex:::.expDataToRanks(test[[1]]$expr)
test1_lab <- test[[1]]$lab
test2_rank <- exprsex:::.expDataToRanks(test[[2]]$expr)
test2_lab <- test[[2]]$lab

ss1_rank <- exprsex:::.expDataToRanks(ss[[1]]$expr)
ss1_lab <- ss[[1]]$lab
ss2_rank <- exprsex:::.expDataToRanks(ss[[2]]$expr)
ss2_lab <- ss[[2]]$lab


sl_genes <- read_clean_sl_genes(organism)
load(file=sprintf("data/04_fits/fit_%s_0304.RData", organism))

# ---------------------------#



# ok now we want the counts of missingness by column
sl_overlap <-intersect(c(sl_genes$f, sl_genes$m), rownames(train_rank))
train_sl <- train_rank[sl_overlap,] # there are 31 total
nas.by.col <- apply(train_sl, 2, function(x) sum(is.na(x)))
summary(nas.by.col)
table(nas.by.col<=1) # 3242 <-- we can try these!
load(file=sprintf("data/04_fits/fit_%s_0304.RData", organism))
preds.1na <- predSexLab(fit, as.matrix(train_sl[,c(nas.by.col<=1)]))
calcAcc(train_lab[c(nas.by.col<=1)], preds.1na) # 96.5%

# missingness by platform??

dat <- train_sl[,c(nas.by.col<=1)]
dat_labels <- train_lab[c(nas.by.col<=1)]

set.seed(116)
removeGenes <-function(nGenes, dat, dat_labels) {
  to.remove <- sl_overlap[sample(1:length(sl_overlap), nGenes, replace=FALSE)]
  dat[to.remove,] <- NA
  acc <- calcAcc(dat_labels, predSexLab(fit, as.matrix(dat)))
  return(data.frame("ngenes"=nGenes, "acc"=acc))
}

## // TODO -- do this with the testing data
acc_remove <- do.call(rbind, lapply(1:50, function(i)
                     do.call(rbind, 
                      lapply(1:length(sl_overlap), 
                             function(x) removeGenes(x, dat, dat_labels)))))

# change NAs to zero
acc_remove[is.na(acc_remove)] <- 0

acc_remove %>% write_csv("data/05_performance/human_missing_train_0204.csv")
# PLOT
ggplot(acc_remove, aes(y=acc, x=factor(ngenes)))+geom_boxplot()+ylab("accuracy (50 runs)")+
  ylim(0,1.0)+xlab("number of genes removed")+
  ggtitle("Sensitivity to gene dropout")
ggsave(sprintf("figures/dropout_%s_train.png", organism), dpi="print", width=6, height=4)

# look at the test data
test_rank <- cbind(test1_rank, test2_rank)
test_lab <- c(test1_lab, test2_lab)
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
acc_remove_test %>% write_csv("data/05_performance/human_missing_test_0204.csv")

# PLOT
ggplot(acc_remove_test, aes(y=acc, x=factor(ngenes)))+geom_boxplot()+ylab("accuracy (50 runs)")+
  ylim(0,1.0)+xlab("number of genes removed")+
  ggtitle("Sensitivity to gene dropout")
ggsave(sprintf("figures/dropout_%s_test.png", organism), dpi="print", width=6, height=4)


# // TODO: 
# - accuracy vs number of genes missing for the test data?
# - gene importance

