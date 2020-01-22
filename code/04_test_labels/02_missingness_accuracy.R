# the goal of this is to assess how much missingness is ok

# look at the the training data from the test platforms


require('exprsex')
require('tidyverse')
require('data.table')


# <--- READ IN THE INPUT ---> #
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

# todo: hard-coded
test.gpls <- c("GPL17586", "GPL10558", "GPL17077", "GPL6480","GPL6102", "GPL4133", "GPL6884", "GPL14550")

gse_gpl2 <- gse_info %>% 
  filter(gpl %in% test.gpls) %>% 
  select(gse, gpl) %>% unique()

keep.dat <- c(train_lab$gse %in% gse_gpl2$gse)
train_rank2 <- train_rank[,keep.dat]
train_sex_lab2 <- train_sex_lab[keep.dat]
save(train_rank2, train_sex_lab2, file="data/03_silver_std/human/03_out_mat/train_plat_filt.RData")
fit2 <- exprsex::trainSexLab(as.matrix(train_rank2),
                            train_sex_lab2,
                            female_genes=fgenes_df2$id,
                            male_genes=mgenes_df2$id)
pred_train2 <- predSexLab(fit2, as.matrix(train_rank2))

confMat(train_sex_lab2, pred_train2)
calcAcc(train_sex_lab2, pred_train2) # 95.1%

# ok now we want the counts of missingness by column
train_sl <- train_rank2[sl_genes,] # there are 30 total
nas.by.col <- apply(train_sl, 2, function(x) sum(is.na(x)))
summary(nas.by.col)
table(nas.by.col==1) # 275 <-- we can try these!
preds.1na <- predSexLab(fit2, as.matrix(train_sl[,c(nas.by.col==1)]))
calcAcc(train_sex_lab2[c(nas.by.col==1)], preds.1na) # 97.1%

dat <- train_sl[,c(nas.by.col==1)]
test_labels <- train_sex_lab2[c(nas.by.col==1)]

set.seed(116)
removeGenes <-function(nGenes, dat, test_labels) {
  to.remove <- sl_genes[sample(1:length(sl_genes), nGenes, replace=FALSE)]
  dat[to.remove,] <- NA
  acc <- calcAcc(test_labels, predSexLab(fit2, as.matrix(dat)))
  return(data.frame("ngenes"=nGenes, "acc"=acc))
}

## // TODO -- do this with the testing data
acc_remove <- do.call(rbind, lapply(1:50, function(i)
                     do.call(rbind, 
                      lapply(1:30, function(x) removeGenes(x, dat, test_labels)))))

# change NAs to zero
acc_remove[is.na(acc_remove)] <- 0

# PLOT
ggplot(acc_remove, aes(y=acc, x=factor(ngenes)))+geom_boxplot()+ylab("accuracy (50 runs)")+
  ylim(0,1.0)+xlab("number of genes removed")+
  ggtitle("Sensitivity to gene dropout")



# // TODO: accuracy vs number of genes missing for the test data?


