# train platform-specific models

require('exprsex')
require('tidyverse')
require('data.table')
source("code/utils/assessment_utils.R")
options(stringsAsFactors = FALSE)

organism <- "human"
run_v <- "full"
# ---- prepare data --- #
common.genes <- read.csv(sprintf("data/consensus_genes_%s.csv", organism))$consensus.genes
if (organism=="human"){
  full.genes <- read.csv(sprintf("data/consensus_genes_%s_full.csv", organism))$consensus.genes
  genes <- full.genes
} else {
  genes <- common.genes
}


train_rank <- read_process_mat(sprintf("data/03_silver_std/%s/03_out_mat/train_%s_rank.csv", organism, run_v), genes)
train_lab <- read_process_lab(sprintf("data/03_silver_std/%s/03_out_mat/train_%s_lab.csv", organism, run_v))

list.test_c <- lapply(0:2, function(i) readInTest(organism, i, "common", common.genes)) 
test_rank_c <- do.call(cbind, lapply(list.test_c, function(x) x$rank))
test_lab_c <- unlist(sapply(list.test_c, function(x) x$lab))

list.test_f <- lapply(0:1, function(i) readInTest(organism, i, "full", full.genes)) 
test_rank_f <- do.call(cbind, lapply(list.test_f, function(x) x$rank))
test_lab_f <- unlist(sapply(list.test_f, function(x) x$lab))

sl_genes <- read_clean_sl_genes(organism, "full")

train_gses <- getGSEsDotted(names(train_lab))
test_c_gses <- getGSEsDotted(names(test_lab_c))
test_f_gses <- getGSEsDotted(names(test_lab_f))

ref <- read_csv(sprintf("data/01_sample_lists/silver_std_%s_reform.csv", organism))
gpls_train <- ref %>% filter(gse %in% train_gses) %>% select(gpl) %>% unique() # gpls
gpls_test_f <- ref %>% filter(gse %in% test_f_gses) %>% select(gpl) %>% unique() # gpls
gpls_test_c <- ref %>% filter(gse %in% test_c_gses) %>% select(gpl) %>% unique() # gpls

# can we save it?

# For each platform:
#   - how many genes are missing?
#   - train a fit
#   - testing and training accuracy

platModels <- function(my.gpl, test_rank, test_lab){
  filt_gpl <- (ref %>% filter(gpl==my.gpl))$gse
  train_keep <- (train_gses %in% filt_gpl)
  train_rank_p <- train_rank[,train_keep]
  train_lab_p <- train_lab[train_keep]
  
  test_gses <- getGSEsDotted(names(test_lab))
  
  test_keep <- (test_gses %in% filt_gpl)
  test_rank_p <- test_rank[,test_keep]
  test_lab_p <- test_lab[test_keep]
  
  fit_p <- exprsex::trainSexLab(train_rank_p, train_lab_p, female_genes = sl_genes$f, male_genes=sl_genes$m)
  pred_train_p <- exprsex::predSexLab(fit_p, train_rank_p)
  pred_test_p <- exprsex::predSexLab(fit_p, test_rank_p)
  
  train_num <- length(unique(sapply(colnames(train_rank_p), function(x) strsplit(x, "\\.")[[1]][[1]])))
  test_num <-  length(unique(sapply(colnames(test_rank_p), function(x) strsplit(x, "\\.")[[1]][[1]])))
  
  train_acc_p <- calcAcc(train_lab_p, pred_train_p)
  test_acc_p <- calcAcc(test_lab_p, pred_test_p)
  
  # instead -- we want the number -CORRECT- by study
  
  fit <- fit_p
  save(fit, file=sprintf("data/04_fits/fit_%s.RData", my.gpl))
  return(data.frame("gpl"=my.gpl, 
                    "studies_train"=train_num, "train_samples"=ncol(train_rank_p),
                    "studies_test"=test_num, "test_samples"=ncol(test_rank_p), "train_acc"=train_acc_p, 
                    "test_acc"=test_acc_p))
}

plat_c <- lapply(gpls_test_c$gpl, platModels(my.gpl, test_rank_c, test_lab_c))
plat_c_df <- do.call(rbind, plat_c)
plat_f <- lapply(gpls_test_f$gpl, function(my.gpl) platModels(my.gpl, test_rank_f, test_lab_f))
plat_f_df <- do.call(rbind, plat_f)
write_csv(rbind(plat_c_df, plat_f_df), "data/05_performance/human_pred_acc_plat_fit.csv")

test_c_o <- read_csv("data/05_performance/human_pred_acc_by_gpl_test_common.csv")
c_df <- right_join(test_c_o %>% select(gpl, acc), plat_c_df %>% 
                     select(gpl, test_acc))  %>% rename(general_fit=acc, plat_fit=test_acc)

test_f_o <- read_csv("data/05_performance/human_pred_acc_by_gpl_test_full.csv")
f_df <- right_join(test_f_o %>% select(gpl, acc), plat_f_df %>% select(gpl,test_acc)) %>% rename(general_fit=acc, plat_fit=test_acc)
