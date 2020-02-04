
# ok so to get a CI for a particular threshold - we can get a p-val

# get the distribution of f-scores for f and for m
# --> construct a 95% CI


require('exprsex')
require('tidyverse')


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
