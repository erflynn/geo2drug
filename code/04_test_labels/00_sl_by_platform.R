# the goal of this script is to set up labels using the train/test data

require('exprsex')
require('tidyverse')
require('data.table')

options(stringsAsFactors = FALSE)

train_rank <- data.table::fread("data/03_silver_std/human/03_out_mat/train_rank.csv", data.table=FALSE)
train_expr <- data.table::fread("data/03_silver_std/human/03_out_mat/train_expr.csv", data.table=FALSE)
train_lab <- read.csv("data/03_silver_std/human/03_out_mat/train_lab.csv")
fgenes_df <- read.csv("data/03_silver_std/human/04_meta_res/fgenes.csv")
mgenes_df <- read.csv("data/03_silver_std/human/04_meta_res/mgenes.csv")

consensus.genes <- read.csv("data/consensus_genes_human.csv")
consensus <- consensus.genes %>% mutate(genes=as.character(consensus.genes))
gse_info <- read_csv("data/01_sample_lists/human_gse_info.csv")



rownames(train_rank) <- consensus$genes
rownames(train_expr) <- consensus$genes

train_lab <- train_lab %>% mutate("sample"=paste(gse, gsm, sep="."))

train_sex_lab <- train_lab$sex
names(train_sex_lab) <- train_lab$sample


readInTest <- function(idx){
  # read in the rest of the test data #
  test_rank <- data.table::fread(sprintf("data/03_silver_std/human/03_out_mat/test_rank_%s.csv", idx), data.table=FALSE)
  test_lab <- read.csv(sprintf("data/03_silver_std/human/03_out_mat/test_lab_%s.csv", idx)) # HUH. these do not match :(
  print(dim(test_rank))
  print(dim(test_lab))
  rownames(test_rank) <- consensus$genes
  test_lab <- test_lab %>% mutate("sample"=paste(gse, gsm, sep=".")) # replace GPL
  test_sex_lab <- test_lab$sex
  names(test_sex_lab) <- test_lab$sample
  test_rank <- test_rank[,names(test_sex_lab)]
  return(list("rank"=test_rank, "lab"=test_sex_lab))
}

# // TODO look into why test_lab is shorter and fix this
list.test <- lapply(0:4, readInTest)

fgenes_df2 <- fgenes_df %>% 
  mutate(id=as.character(entrezgene_id)) %>% 
  select(id) %>%
  unique() 

mgenes_df2 <- mgenes_df %>% 
  mutate(id=as.character(entrezgene_id)) %>% 
  select(id) %>%
  unique() 

confMat <- function(actual, pred){
  table(data.frame(cbind(actual, pred)))
}


calcAcc <- function(actual, pred){
  sum(pred==actual)/length(actual)
}

fit <- exprsex::trainSexLab(as.matrix(train_rank),
                     train_sex_lab,
                     female_genes=fgenes_df2$id,
                     male_genes=mgenes_df2$id)
save(fit, file="data/fit_0120.RData")
pred_train <- predSexLab(fit, as.matrix(train_rank))
pred_test0 <- predSexLab(fit, as.matrix(list.test[[0]]$rank))

confMat(train_sex_lab, pred_train)
confMat(list.test[[0]]$lab, pred_test0)

calcAcc(train_sex_lab, pred_train)
calcAcc(list.test[[0]]$lab, pred_test0)

pred_train_df <- predSexLab(fit, as.matrix(train_rank), scores=TRUE)
pred_train2 <- data.frame(t(pred_train_df))
pred_train2$true_sex <- train_sex_lab
ggplot(pred_train2, aes(x=score_f, y=score_m))+geom_point(aes(color=true_sex, alpha=0.3))+
  geom_abline(slope=1, intercept=fit$threshold) # TODO - add platform here
# this does not look good
#score_m  = fit$threshold + score_f

## BACKGROUND DATA ON PLATFORM ##
train_gse <- gse_info %>% 
  filter(gse %in% unique(str_replace_all(train_lab$gse, "\\.", "-"))) 

pred_train2$gsm <- sapply(rownames(pred_train2 ), function(x) {l <- strsplit(x, "\\.")[[1]]; l[[length(l)]]})
pred_train2$gse <- sapply(rownames(pred_train2 ), function(x) {
  z <- strsplit(x, "\\.")[[1]]; ifelse(length(z)==2, z[[1]], paste(z[[1]], z[[2]], sep="-"))})

pred_train_gpl <- left_join(pred_train2, train_gse %>% select(gse, gpl))

## we see that much of this *IS* separable using the score metric, 
##    but there is large platform-to-platform variability
ggplot(pred_train_gpl, aes(x=score_f, y=score_m))+geom_point(aes(color=true_sex, alpha=0.3))+
  geom_abline(slope=1, intercept=fit$threshold)+facet_wrap(vars(gpl))



# LOOK AT THE TEST DATA
acc_test <- lapply(list.test, function(x) {preds <- predSexLab(fit, as.matrix(x$rank));
calcAcc(x$lab, preds)})

pred_test_df <- predSexLab(fit, as.matrix(list.test[[4]]$rank), scores=TRUE)
pred_test <- data.frame(t(pred_test_df))
pred_test$true_sex <- list.test[[4]]$lab
ggplot(pred_test, aes(x=score_f, y=score_m))+
  geom_point(aes(color=true_sex), alpha=0.5)+ 
  geom_abline(slope=1, intercept=fit$threshold)


getGSEsDotted <- function(list.names){
  # this deals with the fact that many inputted are dotted
  ##samples <- sapply(list.names, function(x) {l <- strsplit(x, "\\.")[[1]]; l[[length(l)]]})
  studies <- sapply(list.names, function(x) {
    z <- strsplit(x, "\\.")[[1]]; ifelse(length(z)==2, z[[1]], paste(z[[1]], z[[2]], sep="-"))})
  return(studies)
}


byPlatAcc <- function(fit, df, lab, my.gpl){
  gses <- unique(sapply((gse_info %>% filter(gpl == my.gpl))$gse, function(x) strsplit(x, "\\.")[[1]][[1]]))
  to.keep <- sapply(colnames(df), function(x) ifelse(strsplit(x, "\\.")[[1]][[1]] %in% gses, TRUE, FALSE))
  # get the overlap in gses2
  df.gpl <- df[,to.keep]
  gses2 <- unique(lapply(colnames(df.gpl), function(x) strsplit(x, "\\.")[[1]][[1]]))
  preds <- predSexLab(fit, as.matrix(df.gpl))
  acc <- calcAcc(preds, lab[to.keep])
  return(data.frame("gpl"=my.gpl, "acc"=acc, "num_studies"=length(gses2), "num_samples"=length(preds)))
}

all_plats_res <- function(fit, df, lab){
  gses <- unique(getGSEsDotted(colnames(df)))
  gpl_lst <- unique((gse_info %>% filter(gse %in% gses))$gpl)
  do.call(rbind, lapply(gpl_lst, function(my.gpl) byPlatAcc(fit, df, lab, my.gpl)))
}
train_acc <- all_plats_res(fit, train_rank, train_sex_lab) %>% arrange(acc)

# ---- we need to group the test data by GPL ---- #
# there are only eight --> not so bad

# get a list of test gpls
test.gses <- unique(unlist(lapply(list.test, function(x) getGSEsDotted(names(x$lab)))))
test.info <-  gse_info %>% filter(gse %in% test.gses) 
test.gpls <- (test.info %>% select(gpl) %>% unique())$gpl


by_plat_test <- function(list.dat, my.gpl){
  gse_gpl <- test.info %>% select(gse, gpl) %>% filter(gpl==my.gpl) %>% select(gse) %>% unique()
  list.dat2 <- lapply(list.dat, function(dat){
    lab <- dat$lab; 
    df <- dat$rank;
    gses.dat <- getGSEsDotted(names(lab))
    to.keep <- c(gses.dat %in% gse_gpl$gse)
    lab2 <- lab[to.keep]
    df2 <- df[,to.keep]
    return(list("lab"=lab2, "rank"=df2))
  })
  rank2 <- do.call(cbind, lapply(list.dat2, function(x) x$rank))
  lab2 <- unlist(lapply(list.dat2, function(x) x$lab))
  return(list("rank"=rank2, "lab"=lab2))
}

lapply(test.gpls, function(my.gpl){
  test_gpl <- by_plat_test(list.test, my.gpl)
  save(test_gpl, file=sprintf("data/03_silver_std/human/03_out_mat/test_%s.RData", tolower(my.gpl)))
  
})

# by platform-accuracy in the test data
test_plat_acc <- function(my.gpl){
  print(my.gpl)
  miceadds::load.Rdata(sprintf("data/03_silver_std/human/03_out_mat/test_%s.RData", tolower(my.gpl)), "test_gpl")
  preds_test_df <- predSexLab(fit, as.matrix(test_gpl$rank), scores=TRUE)
  num_samples <- length(test_gpl$lab)
  num_studies <- length(unique(sapply(names(test_gpl$lab), function(x) strsplit(x, "\\.")[[1]][[1]])))
  
  pred_test <- data.frame(t(pred_test_df))
  pred_test$true_sex <- test_gpl$lab
  ggplot(pred_test, aes(x=score_f, y=score_m))+
    geom_point(aes(color=true_sex), alpha=0.5)+ 
    geom_abline(slope=1, intercept=fit$threshold)+
    ggtitle(sprintf("%s (test)", my.gpl))
  ggsave(file=sprintf("test_score_%s.png", tolower(my.gpl)), width=5, height=4)
  acc_gpl <- calcAcc(pred_test$sex, test_gpl$lab)
  
  return(data.frame("gpl"=my.gpl, "acc"=acc_gpl, "num_studies"=num_studies, "num_samples"=num_samples))
}

# plot for test data

test_acc <- lapply(test.gpls, test_plat_acc)
test_acc_df <- do.call(rbind, test_acc)

#### TODO - save the results

train_acc %>% filter(gpl %in% test.gpls) %>% summarize(ntot=sum(num_studies), nsamp=sum(num_samples))
# 13 studies, 808 samples, 8 platforms

train_acc %>% filter(!gpl %in% test.gpls) %>% arrange(desc(num_studies)) 
# 3 of these have 2 -- so we *COULD* test on better
# but what about the rest??


## // TODO - this does not work because of NAs?? or b/c we need to adjust??
#fit_expr <- exprsex::trainSexLab(as.matrix(train_expr),
#                            train_sex_lab,
#                            female_genes=fgenes_df2$id,
#                            male_genes=mgenes_df2$id)
