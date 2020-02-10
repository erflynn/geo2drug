# the goal of this script is to set up labels using the train/test data

require('exprsex')
require('tidyverse')
require('data.table')
source("code/utils/assessment_utils.R")
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]

# ---- load the data ---- #
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


# ---- train and save a model ---- #
fit <- exprsex::trainSexLab(train_rank,
                            train_lab,
                            female_genes=sl_genes$f,
                            male_genes=sl_genes$m)
suffix <- "0304"
save(fit, file=sprintf("data/04_fits/fit_%s_%s.RData", organism, suffix))
pred_train <- predSexLab(fit, train_rank)
pred_test <- predSexLab(fit, test1_rank)
pred_test2 <- predSexLab(fit, test2_rank)

pred_ss <- predSexLab(fit, ss1_rank)
pred_ss2 <- predSexLab(fit, ss2_rank)

confMat(train_lab, pred_train)
calcAcc(train_lab, pred_train) # 94.7% training accuracy 

confMat(pred_test, test1_lab)
calcAcc(pred_test, test1_lab ) # 94.8% test accuracy 

confMat(pred_test2, test2_lab)
calcAcc(pred_test2, test2_lab ) # 74.5% test accuracy -- lower 

confMat(c(pred_test, pred_test2), c(test1_lab, test2_lab))
calcAcc(c(pred_test, pred_test2), c(test1_lab, test2_lab)) # 88.6%


confMat(pred_ss, ss1_lab)
calcAcc(pred_ss, ss1_lab ) # 82.2% accuracy 
confMat(pred_ss2, ss2_lab)
calcAcc(pred_ss2, ss2_lab ) # 94.1% accuracy

confMat(c(pred_ss, pred_ss2), c(ss1_lab, ss2_lab))
calcAcc(c(pred_ss, pred_ss2), c(ss1_lab, ss2_lab)) # 86.0% overall
# <--- I think we have a lot of platform artifacts going on :(

# ok so I'm not sure *IF* the single-sex studies have lower accuracy 
#  -OR- if this is a platform-specific artifact
# --> need to plot platform accuracy broken down by category
# first step:
#  - calculate by study accuracy

# --- grab gpl info --- #
ss_ref <- read_csv("data/01_sample_lists/human_single_sex_studies.csv")
ss_ref2 <- ss_ref %>% sepReformatGPL()
silver_std <- read_csv(sprintf("data/01_sample_lists/silver_std_%s_reform.csv", organism))
train_pred_df <- pred_out_df(fit, train_rank, train_lab, organism, silver_std)
train_rm <- train_pred_df %>% group_by(gse) %>% summarize(num_f=sum(true_sex==1), num_m=sum(true_sex==0)) %>%
  filter(num_f==0 | num_m==0)
train_pred_df2 <- train_pred_df %>% filter(!gse %in% train_rm$gse)
calcAcc(train_pred_df2$pred_sex, train_pred_df2$true_sex)
train_stud <- acc_by_study(train_pred_df2)

ss_pred_df <- pred_out_df(fit, cbind(ss1_rank, ss2_rank), 
                          c(ss1_lab, ss2_lab), organism, ss_ref2)
ss_stud <- acc_by_study(ss_pred_df)

test_pred_df <- pred_out_df(fit, cbind(test1_rank, test2_rank), 
                          c(test1_lab, test2_lab), organism, silver_std)

test_rm <- test_pred_df %>% group_by(gse) %>% summarize(num_f=sum(true_sex==1), num_m=sum(true_sex==0)) %>%
  filter(num_f==0 | num_m==0)
test_pred_df2 <- test_pred_df %>% filter(!gse %in% test_rm$gse)
calcAcc(test_pred_df2$pred_sex, test_pred_df2$true_sex)
test_stud <- acc_by_study(test_pred_df2)

# --- plot the accuracy --- #
test_stud$dataset <- "test"
train_stud$dataset <- "train"
ss_stud <- ss_stud %>% filter(gpl %in% c(unique(train_stud$gpl)))
ss_stud$dataset <- "single_sex"
ds_stud_acc <- test_stud %>% bind_rows(train_stud) %>% bind_rows(ss_stud)
ds_stud_acc$dataset <- factor(ds_stud_acc$dataset, levels=c("train", "test", "single_sex"))
train_gpl_acc_order <- train_stud %>% group_by(gpl) %>% summarize(acc=mean(acc)) %>% arrange(desc(acc))
ds_stud_acc$gpl <- factor(ds_stud_acc$gpl, levels=train_gpl_acc_order$gpl)
ggplot(ds_stud_acc, aes(x=gpl, y=acc))+
  geom_boxplot(aes(col=dataset))+geom_point(aes(col=dataset), position=position_dodge(0.8))+
  theme(axis.text.x = element_text(angle = 90))+ylab("accuracy")+xlab("")
ggsave(file="figures/plat_acc_w_ss.png", dpi="print", width=6, height=3.5)

# ---- visualize the SS data by platform ---- #
summarizeAcc2 <- function(pred_df){
  pred_acc <- pred_df %>% group_by(gpl) %>% 
    summarize(num_studies=length(unique(gse)), 
              true_pred=sum(true_sex==pred_sex, na.rm=TRUE), num_samples=length(true_sex)) %>%
    mutate(acc=true_pred/num_samples)
  pred_overall <- pred_df %>% 
    summarize(num_studies=length(unique(gse)), true_pred=sum(true_sex==pred_sex, na.rm=TRUE), 
              num_samples=length(true_sex)) %>%
    mutate(acc=true_pred/num_samples) %>%
    mutate(gpl="overall") %>% select(gpl, everything())
  pred_acc2 <- rbind(pred_acc, pred_overall)
  return(pred_acc2)
}

train_gpl <- summarizeAcc2(train_pred_df2)
test_gpl <- summarizeAcc2(test_pred_df2)
ss_gpl <- summarizeAcc2(ss_pred_df)
ggplot(train_pred_df2 %>% semi_join(train_stud %>% filter(acc < 0.9), by="gse"), aes(x=score_f, y=score_m))+geom_point(aes(color=true_sex, alpha=0.3))+
  geom_abline(slope=1, intercept=fit$threshold)+facet_wrap(vars(gpl,gse))
ggplot(test_pred_df2 %>% semi_join(test_stud %>% filter(acc < 0.9), by="gse"), aes(x=score_f, y=score_m))+geom_point(aes(color=true_sex, alpha=0.3))+
  geom_abline(slope=1, intercept=fit$threshold)+facet_wrap(vars(gpl, gse))
ggsave("figures/test_poor_perf.png")

ggplot(ss_pred_df %>% semi_join(ss_stud %>% filter(acc < 0.9), by="gse"), aes(x=score_f, y=score_m))+geom_point(aes(color=true_sex, alpha=0.3))+
  geom_abline(slope=1, intercept=fit$threshold)+facet_wrap(vars(gpl, gse))
ggsave("figures/ss_poor_perf.png")
# eda on the terrible studies - what is going here?


# ---  test the model --- #
#  - compute the full test accuracy
#  - compute the by-platform test accuracy



# if (organism == "human"){
#   # // TODO look into why test_lab is shorter and fix this
#   list.test_c <- lapply(0:2, function(i) readInTest(organism, i, "common", common.genes)) 
#   list.test_f <- lapply(0:1, function(i) readInTest(organism, i, "full", full.genes))
#   list.test_r <- lapply(0:2, function(i) readInTest(organism, i, "rare", full.genes))
#   test_c_pred <- do.call(rbind, lapply(1:3, function(i) 
#                            pred_out_df(fit, list.test_c[[i]]$rank, list.test_c[[i]]$lab, organism)))
#   
#   test_f_pred <- do.call(rbind, lapply(1:2, function(i) 
#                            pred_out_df(fit, list.test_f[[i]]$rank, list.test_f[[i]]$lab, organism)))
#   
#   test_r_pred <- do.call(rbind, lapply(1:3, function(i) 
#                            pred_out_df(fit, list.test_r[[i]]$rank, list.test_r[[i]]$lab, organism)))
#   
#   train_common <- read_process_mat(sprintf("data/03_silver_std/%s/03_out_mat/train_common_rank.csv", organism), common.genes)
#   train_c_lab <- read_process_lab(sprintf("data/03_silver_std/%s/03_out_mat/train_common_lab.csv", organism))
#   train_c_pred <- pred_out_df(fit, train_common, train_c_lab, organism)
#   
#   summarizeAcc(train_c_pred, "train_common", organism)
#   summarizeAcc(test_c_pred, "test_common", organism)
#   summarizeAcc(rbind(test_c_pred, test_f_pred), "test_full", organism)
#   summarizeAcc(test_f_pred, "test_full_sm", organism)
#   summarizeAcc(test_r_pred, "test_rare", organism)
# } else{
#   test_c <-  readInTest(organism, 0, "common", common.genes)
#   test_c_pred <- pred_out_df(fit, test_c$rank, test_c$lab, organism)
#   summarizeAcc(test_c_pred, "test", organism)
#   
#   if (organism == "mouse"){
#     test_r <-  readInTest(organism, 0, "rare", common.genes)
#     test_r_pred <- pred_out_df(fit, test_r$rank, test_r$lab, organism)
#     summarizeAcc(test_r_pred, "test_rare", organism)
#     
#   }
# } 
# 




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
