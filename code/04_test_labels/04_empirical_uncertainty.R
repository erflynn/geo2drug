
# ok so to get a probability for each prediction / measurement of uncertainty

# get the distribution of f-scores for f and for m
# --> get a proba score by fitting a logit / sigmoid function
# HOW DID I CREATE ROC CURVES?


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
suffix <- 0304
load(fit, file=sprintf("data/04_fits/fit_%s_%s.RData", organism, suffix))

# --- grab gpl info --- #
ss_ref <- read_csv("data/01_sample_lists/human_single_sex_studies.csv")
ss_ref2 <- ss_ref %>% sepReformatGPL()
silver_std <- read_csv(sprintf("data/01_sample_lists/silver_std_%s_reform.csv", organism))
train_pred_df <- pred_out_df(fit, train_rank, train_lab, organism, silver_std)
train_rm <- train_pred_df %>% group_by(gse) %>% summarize(num_f=sum(true_sex==1), num_m=sum(true_sex==0)) %>%
  filter(num_f==0 | num_m==0)
train_pred_df2 <- train_pred_df %>% filter(!gse %in% train_rm$gse)

ss_pred_df <- pred_out_df(fit, cbind(ss1_rank, ss2_rank), 
                          c(ss1_lab, ss2_lab), organism, ss_ref2)

test_pred_df <- pred_out_df(fit, cbind(test1_rank, test2_rank), 
                            c(test1_lab, test2_lab), organism, silver_std)

test_rm <- test_pred_df %>% group_by(gse) %>% summarize(num_f=sum(true_sex==1), num_m=sum(true_sex==0)) %>%
  filter(num_f==0 | num_m==0)
test_pred_df2 <- test_pred_df %>% filter(!gse %in% test_rm$gse)

# how did a I do a roc plot before?
head(test_pred_df2)

test_pred_df3 <- test_pred_df2 %>% mutate(diff_score=(score_m-score_f))
scores <- test_pred_df3$diff_score
roc.r <- roc(test_pred_df3$true_sex,scores, plot=TRUE)
roc.r$auc

ss_pred_df3 <- ss_pred_df %>% mutate(diff_score=(score_m-score_f))
roc.r <- roc(ss_pred_df3$true_sex,ss_pred_df3$diff_score, plot=TRUE)
roc.r$auc

require('qrnn')
train.scores <- train_pred_df2$score_m-train_pred_df2$score_f
sigmoid(train.scores)
plot(train.scores, sigmoid(train.scores), col=ifelse(train_pred_df2$true_sex==0, "red", "blue"))

require('sigmoid')
train.scores <- train_pred_df2$score_m-train_pred_df2$score_f
train_pred_df2 <- train_pred_df2 %>% mutate(diff_score=(score_m-score_f), true_sex=as.factor(true_sex))

plot(train.scores, SoftMax(train.scores), col=ifelse(train_pred_df2$true_sex==0, "red", "blue"))

lm.fit <- glm(true_sex ~ diff_score, data=train_pred_df2, family="binomial")
train.pred <- predict(lm.fit, train_pred_df2, type="response")
train_pred_df2$train.pred <- train.pred
ggplot(train_pred_df2, aes(x=true_sex, y=train.pred))+geom_violin(aes(color=true_sex))
roc(train_pred_df2$true_sex, train_pred_df2$train.pred) # 0.988

sapply(seq(0,1, 0.05), function(x){
  train_pred_df2$pred_class <- ifelse(train_pred_df2$train.pred > x, 1, 0)
  return(calcAcc(train_pred_df2$true_sex,train_pred_df2$pred_class))
})
# in between threshold

roc.df <- roc_(seq(0,1, 0.05), train_pred_df2$true_sex, train_pred_df2$train.pred, ret="coords")
train_pred_df2$pred_class <- ifelse(train_pred_df2$train.pred > 0.5, 1, 0)

calcAcc(train_pred_df2$true_sex,train_pred_df2$pred_class)

test_pred_df3$pred <- predict(lm.fit, test_pred_df3, type="response")
ggplot(test_pred_df3, aes(x=factor(true_sex), y=pred))+geom_violin(aes(color=factor(true_sex)))
roc(test_pred_df3$true_sex, test_pred_df3$pred, plot=TRUE) # 0.976
test_pred_df3$pred_class <- ifelse(test_pred_df3$pred > 0.5, 1, 0)
sapply(seq(0,1, 0.05), function(x){
  test_pred_df3$pred_class <- ifelse(test_pred_df3$pred > x, 1, 0)
  return(calcAcc(test_pred_df3$true_sex,test_pred_df3$pred_class))
})

test_pred_df4 <- test_pred_df3 %>% filter(pred > 0.9)
confMat(test_pred_df3$true_sex,test_pred_df3$pred_class)
calcAcc(test_pred_df4$true_sex,test_pred_df4$pred_class) 

plot(test_pred_df3$diff_score, test_pred_df3$pred, col=ifelse(test_pred_df3$true_sex==0, "red", "blue"))

sapply(c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2), function(x){
test_pred_df4 <- test_pred_df3 %>% filter(pred > 1-x | pred < x)
calcAcc(test_pred_df4$true_sex,test_pred_df4$pred_class) 
})
# accuracy is that or greater. we are fine

###### DEPRECATED ####### 

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
test_rank <- do.call(cbind, lapply(list.test_c, function(x) x$rank))
test_lab <- unlist(sapply(list.test_c, function(x) x$lab))

load(file=sprintf("data/04_fits/fit_%s_%s.RData", organism, run_v))
pred <- predSexLab(fit, as.matrix(train_rank), scores=TRUE)
pred <- pred %>% rename(pred_sex=sex) 
pred$actual_sex <- train_lab
pred <- pred %>% mutate(diff_score=score_m-score_f)
calcAcc(pred$pred_sex, pred$actual_sex) # 0.952
mean_sd <- pred %>% select(-pred_sex) %>% 
  pivot_longer(cols=c(score_m, score_f, diff_score), names_to="score") %>%
  group_by(actual_sex, score) %>%
  summarize(mean=mean(value), sd=sd(value))
#%>%
#  mutate(ci_95_u=mean+1.96*sd, ci_95_l=mean-1.96*sd)
upper_bound_f <- mean_sd %>% filter(actual_sex==0 & score=="diff_score") %>% mutate(val=mean+1.64*sd)
lower_bound_m <- mean_sd %>% filter(actual_sex==1 & score=="diff_score") %>% mutate(val=mean-1.64*sd)

# we want 1.64*sd here not 1.96 --> 95% 1-sided CI (45%s)

pred2 <- pred %>% mutate(pred_sex2=ifelse(
  (actual_sex == 0 & (diff_score > upper_bound_f$val)) |
                     (actual_sex == 1 & (diff_score < lower_bound_m$val)), NA, pred_sex))
calcAcc((pred2 %>% filter(!is.na(pred_sex2)))$pred_sex,
        (pred2 %>% filter(!is.na(pred_sex2)))$actual_sex )
# 98.4% accuracy


### look at the rest of the

#### look at the testing data #####
cutoffScores <- function(pred_df){
  pred_df %>%  mutate(pred_sex2=ifelse((actual_sex == 0 & (diff_score > upper_bound_f$val)) |
    (actual_sex == 1 & (diff_score < lower_bound_m$val)), NA, pred_sex))
}

pred_test <- predSexLab(fit, as.matrix(test_rank), scores=TRUE)
pred_test <- pred_test %>% rename(pred_sex=sex) 
pred_test$actual_sex <- test_lab
pred_test <- pred_test %>% mutate(diff_score=score_m-score_f)
pred_test2 <- cutoffScores(pred_test)
calcAcc(pred_test$pred_sex, pred_test$actual_sex) # 0.934
calcAcc((pred_test2 %>% filter(!is.na(pred_sex2)))$pred_sex,
        (pred_test2 %>% filter(!is.na(pred_sex2)))$actual_sex )


test_unc_filt <- function(my.gpl){
  print(my.gpl)
  miceadds::load.Rdata(sprintf("data/03_silver_std/human/03_out_mat/test_%s.RData", tolower(my.gpl)), "test_gpl")
  preds_test_df <- predSexLab(fit, as.matrix(test_gpl$rank), scores=TRUE)
  preds_test_df <- preds_test_df %>% rename(pred_sex=sex) %>% mutate(diff_score=score_m-score_f)
  preds_test_df$actual_sex <- test_gpl$lab

  preds_test_df2 <- cutoffScores(preds_test_df) 
  num_samples <- length(test_gpl$lab)
  num_studies <- length(unique(sapply(names(test_gpl$lab), 
                                      function(x) strsplit(x, "\\.")[[1]][[1]])))
  
  acc_gpl <- calcAcc(preds_test_df2$pred_sex, preds_test_df2$actual_sex)
  pred3 <- preds_test_df2 %>% filter(!is.na(pred_sex2))
  acc_gpl2 <- calcAcc(pred3$pred_sex, pred3$actual_sex)
  frac.rem <- sum(is.na(preds_test_df2$pred_sex2))/nrow(preds_test_df)
  print(sprintf("%s: %s, %s, %s", my.gpl, acc_gpl, acc_gpl2, frac.rem))
  return(data.frame("gpl"=my.gpl, "acc1"=acc_gpl, "acc2"=acc_gpl2, "frac"=frac.rem, "num_studies"=num_studies, "num_samples"=num_samples))
}

test_acc <- lapply(test.gpls, test_unc_filt)

test_acc_df <- do.call(rbind, test_acc)
save(test_acc_df, file="data/03_silver_std/human/test_acc_rem.RData")



options(stringsAsFactors = FALSE)

train_rank <- data.table::fread("data/03_silver_std/human/03_out_mat/train_rank.csv", data.table=FALSE)
train_lab <- read.csv("data/03_silver_std/human/03_out_mat/train_lab.csv")
train_lab <- train_lab %>% mutate(gse=str_replace_all(gse, "\\.", "-"))
fgenes_df <- read.csv("data/03_silver_std/human/04_meta_res/fgenes.csv")
mgenes_df <- read.csv("data/03_silver_std/human/04_meta_res/mgenes.csv")

fgenes_df2 <- fgenes_df %>% 
  mutate(id=as.character(entrezgene_id)) %>% 
  select(id) %>%
  unique() 

mgenes_df2 <- mgenes_df %>% 
  mutate(id=as.character(entrezgene_id)) %>% 
  select(id) %>%
  unique()
sl_genes <- c(fgenes_df2$id, mgenes_df2$id) 

consensus.genes <- read.csv("data/consensus_genes_human.csv")
consensus <- consensus.genes %>% mutate(genes=as.character(consensus.genes))
rownames(train_rank) <- consensus$genes
train_lab <- train_lab %>% mutate("sample"=paste(gse, gsm, sep="."))
train_sex_lab <- train_lab$sex
names(train_sex_lab) <- train_lab$sample

gse_info <- read_csv("data/01_sample_lists/human_gse_info.csv")
