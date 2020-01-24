# the goal of this script is to set up labels using the train/test data

require('exprsex')
require('tidyverse')
require('data.table')
source("code/utils/assessment_utils.R")
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
run_v <- args[2]



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


# ---- train and save a model ---- #
fit <- exprsex::trainSexLab(train_rank,
                            train_lab,
                            female_genes=sl_genes$f,
                            male_genes=sl_genes$m)
save(fit, file=sprintf("data/04_fits/fit_%s_%s.RData", organism, run_v))
pred_train <- predSexLab(fit, train_rank)

confMat(train_lab, pred_train)
calcAcc(train_lab, pred_train) # 93.9% training accuracy for common


# --- grab gpl info --- #
silver_std <- read.csv(sprintf("data/01_sample_lists/gse_for_silver_std_%s.csv", organism), stringsAsFactors = FALSE)
silver_std_s <- silver_std %>% filter(!str_detect(gpl, "\\|")) 
silver_std_m <- silver_std %>% filter(str_detect(gpl, "\\|"))  %>% separate_rows(gpl, sep="\\|") %>%
  mutate(gse=sprintf("%s-%s", gse, gpl)) 
silver_std2 <- data.frame(rbind(silver_std_s, silver_std_m))
silver_std2 %>% write_csv(sprintf("data/01_sample_lists/silver_std_%s_reform.csv", organism))

train_pred_df <- pred_out_df(fit, train_rank, train_lab, organism)
summarizeAcc(train_pred_df, "train", organism)

# ---  test the model --- #
#  - compute the full test accuracy
#  - compute the by-platform test accuracy
# // TODO: - compute the by-study test accuracy --> plot the accuracy by platform across all

if (organism == "human"){
  # // TODO look into why test_lab is shorter and fix this
  list.test_c <- lapply(0:2, function(i) readInTest(organism, i, "common", common.genes)) 
  list.test_f <- lapply(0:1, function(i) readInTest(organism, i, "full", full.genes))
  list.test_r <- lapply(0:2, function(i) readInTest(organism, i, "rare", full.genes))
  test_c_pred <- do.call(rbind, lapply(1:3, function(i) 
                           pred_out_df(fit, list.test_c[[i]]$rank, list.test_c[[i]]$lab, organism)))
  
  test_f_pred <- do.call(rbind, lapply(1:2, function(i) 
                           pred_out_df(fit, list.test_f[[i]]$rank, list.test_f[[i]]$lab, organism)))
  
  test_r_pred <- do.call(rbind, lapply(1:3, function(i) 
                           pred_out_df(fit, list.test_r[[i]]$rank, list.test_r[[i]]$lab, organism)))
  
  train_common <- read_process_mat(sprintf("data/03_silver_std/%s/03_out_mat/train_common_rank.csv", organism), common.genes)
  train_c_lab <- read_process_lab(sprintf("data/03_silver_std/%s/03_out_mat/train_common_lab.csv", organism))
  train_c_pred <- pred_out_df(fit, train_common, train_c_lab, organism)
  
  summarizeAcc(train_c_pred, "train_common", organism)
  summarizeAcc(test_c_pred, "test_common", organism)
  summarizeAcc(rbind(test_c_pred, test_f_pred), "test_full", organism)
  summarizeAcc(test_f_pred, "test_full_sm", organism)
  summarizeAcc(test_r_pred, "test_rare", organism)
} else{
  test_c <-  readInTest(organism, 0, "common", common.genes)
  test_c_pred <- pred_out_df(fit, test_c$rank, test_c$lab, organism)
  summarizeAcc(test_c_pred, "test", organism)
  
  if (organism == "mouse"){
    test_r <-  readInTest(organism, 0, "rare", common.genes)
    test_r_pred <- pred_out_df(fit, test_r$rank, test_r$lab, organism)
    summarizeAcc(test_r_pred, "test_rare", organism)
    
  }
} 





# ------ scratch work ----- #

# # get a list of test gpls
# test.gses <- unique(unlist(lapply(list.test, function(x) getGSEsDotted(names(x$lab)))))
# test.info <-  gse_info %>% filter(gse %in% test.gses) 
# test.gpls <- (test.info %>% select(gpl) %>% unique())$gpl
# 
# 
# by_plat_test <- function(list.dat, my.gpl){
#   gse_gpl <- test.info %>% select(gse, gpl) %>% filter(gpl==my.gpl) %>% select(gse) %>% unique()
#   list.dat2 <- lapply(list.dat, function(dat){
#     lab <- dat$lab; 
#     df <- dat$rank;
#     gses.dat <- getGSEsDotted(names(lab))
#     to.keep <- c(gses.dat %in% gse_gpl$gse)
#     lab2 <- lab[to.keep]
#     df2 <- df[,to.keep]
#     return(list("lab"=lab2, "rank"=df2))
#   })
#   rank2 <- do.call(cbind, lapply(list.dat2, function(x) x$rank))
#   lab2 <- unlist(lapply(list.dat2, function(x) x$lab))
#   return(list("rank"=rank2, "lab"=lab2))
# }
# 
# lapply(test.gpls, function(my.gpl){
#   test_gpl <- by_plat_test(list.test, my.gpl)
#   save(test_gpl, file=sprintf("data/03_silver_std/human/03_out_mat/test_%s.RData", tolower(my.gpl)))
#   
# })
# 
# # by platform-accuracy in the test data
# test_plat_acc <- function(my.gpl){
#   print(my.gpl)
#   miceadds::load.Rdata(sprintf("data/03_silver_std/human/03_out_mat/test_%s.RData", tolower(my.gpl)), "test_gpl")
#   preds_test_df <- predSexLab(fit, as.matrix(test_gpl$rank), scores=TRUE)
#   num_samples <- length(test_gpl$lab)
#   num_studies <- length(unique(sapply(names(test_gpl$lab), function(x) strsplit(x, "\\.")[[1]][[1]])))
#   
#   pred_test <- data.frame(t(pred_test_df))
#   pred_test$true_sex <- test_gpl$lab
#   ggplot(pred_test, aes(x=score_f, y=score_m))+
#     geom_point(aes(color=true_sex), alpha=0.5)+ 
#     geom_abline(slope=1, intercept=fit$threshold)+
#     ggtitle(sprintf("%s (test)", my.gpl))
#   ggsave(file=sprintf("test_score_%s.png", tolower(my.gpl)), width=5, height=4)
#   acc_gpl <- calcAcc(pred_test$sex, test_gpl$lab)
#   
#   return(data.frame("gpl"=my.gpl, "acc"=acc_gpl, "num_studies"=num_studies, "num_samples"=num_samples))
# }
# 
# # plot for test data
# 
# test_acc <- lapply(test.gpls, test_plat_acc)
# test_acc_df <- do.call(rbind, test_acc)
# 
# train_acc %>% filter(gpl %in% test.gpls) %>% summarize(ntot=sum(num_studies), nsamp=sum(num_samples))
# # 13 studies, 808 samples, 8 platforms
# 
# train_acc %>% filter(!gpl %in% test.gpls) %>% arrange(desc(num_studies)) 
# # 3 of these have 2 -- so we *COULD* test on better
# # but what about the rest??
# 
# 
# ## // TODO - this does not work because of NAs?? or b/c we need to adjust??
# #fit_expr <- exprsex::trainSexLab(as.matrix(train_expr),
# #                            train_sex_lab,
# #                            female_genes=fgenes_df2$id,
# #                            male_genes=mgenes_df2$id)
