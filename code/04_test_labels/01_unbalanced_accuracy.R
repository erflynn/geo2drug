source("code/utils/sex_lab_utils.R")
source("code/utils/assessment_utils.R")

require('exprsex')
require('tidyverse')

## load the data

ss <- read_process_ds(organism, "06_single_sex", "single_sex")
train <- read_process_ds(organism, "03_silver_std", "train_common")
train_full <- read_process_ds(organism, "03_silver_std", "train_full")
test <- read_process_ds(organism, "03_silver_std", "testing")
test_full <- read_process_ds(organism, "03_silver_std", "testing_full")

miceadds::load.Rdata(sprintf("../gpl_ref/%s_xy_genes.RData", organism), "xy_genes") 
xy_genes2 <- xy_genes %>% unique() %>% filter(!is.na(gene))
x_genes <- xy_genes2 %>% filter(chr=="X") %>% dplyr::select(gene)
y_genes <- xy_genes2 %>% filter(chr=="Y") %>% dplyr::select(gene)
xy_genes2.2 <- xy_genes2$gene

miceadds::load.Rdata("../gpl_ref/human_gene_map.RData", "gene_map")
toker_list <- c("XIST", "KDM5D", "RPS4Y1")
toker_list2 <- gene_map %>% filter(hgnc_symbol %in% toker_list) %>% dplyr::select(hgnc_symbol, entrezgene_id) %>% unique()
toker_list2.2 <- toker_list2$entrezgene_id
toker_list2f <- (toker_list2 %>% filter(hgnc_symbol=="XIST") )$entrezgene_id
toker_list2m <- (toker_list2 %>% filter(hgnc_symbol %in% c("KDM5D", "RPS4Y1") ))$entrezgene_id


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


downSampleTestDat <- function(fraction_m, gse, test_expr, test_lab){
  test_cols <- filterColsByGSE(colnames(test_expr), c(gse))
  test_exp1 <- test_expr[,test_cols]
  test_rank1 <- test_rank[,test_cols]
  
  test_labels1 <- test_lab[test_cols]
  
  # -- NOW DOWNSAMPLE THE TESTING DATA -- #
  m_samples <- names(test_labels1)[(test_labels1==1)]
  f_samples <- names(test_labels1)[(test_labels1==0)]
  fsamples <- test_exp1[,f_samples]
  msamples <- test_exp1[,m_samples]
  
  nsamples <- length(test_labels1)
  sex_counts <- table(test_labels1)
  count_f <- sum(test_labels1==0)
  count_m <- sum(test_labels1==1)
  actual_frac_m <- count_m/nsamples
  actual_frac_f <- count_f/nsamples
  if (fraction_m == 0){
    downsampled_df <- fsamples
  }
  else if (fraction_m == 1){
    downsampled_df <- msamples
  } 
  # now downsample males or females
  else if (actual_frac_m > fraction_m){
    num_m <- floor(count_f/(1-fraction_m))-count_f
    if (num_m < 2){ num_m <- 2}
    
    mcols <- sample(ncol(msamples), num_m)
    downsampled_df <- cbind(msamples[,mcols], fsamples)
  } 
  else {
    num_f <- floor(count_m/fraction_m)-count_m
    if (num_f < 2){ num_f <- 2}
    
    fcols <- sample(ncol(fsamples), num_f)
    downsampled_df <- cbind(fsamples[,fcols], msamples)
  }
  
  downsampled_labels <- test_lab[colnames(downsampled_df)]
  new.order <- order(names(downsampled_labels))
  downsampled_labels <- downsampled_labels[new.order]  
  downsampled_df <- downsampled_df[,new.order]
  return(list("df"=downsampled_df, "lab"=downsampled_labels))
}

calcAccDown <- function(fraction_m,
                        train_expr, train_rank, train_lab,
                        test_expr, test_rank, test_lab){
  # calculate the accuracy for downsampled data
  rank_lab <- c()
  actual_lab <- c()
  toker_lab <- c()
  massir_lab <- c()
  gses <- unique(getGSEsDotted(colnames(test_expr)))
  # for (j in 1:length(gpl_list)){
  #   gpl <- gpl_list[[j]]
  #   print(gpl)
  #   
  #   train_cols <- filterColsByGSE(colnames(train_expr), train_gpl_to_gse[[gpl]])
  #   train_exp1 <- train_expr[,train_cols]
  #   train_rank1 <- train_rank[,train_cols]
  #   
  #   train_labels1 <- train_lab[train_cols]
  #   gses <- test_gpl_to_gse[[gpl]]
    
    for (i in 1:length(gses)){
      gse <- gses[i]
      ds <- downSampleTestDat(fraction_m, gse, test_expr, test_lab)
      downsampled_df <- ds$df
      downsampled_labels <- ds$lab
      
      if (ncol(downsampled_df)!=1 & !is.null(colnames(downsampled_df))){
        expr2 <- downsampled_df
        keys <- rownames(expr2)
        list.keys <- keys[keys %in% xy_genes2.2]
        
        # massiR labels
        ychr_keys <- keys[keys %in% y_genes$gene]
        massir_sex_labels <- massiRAcc(expr2, ychr_keys, min.samples=3) 
        if (is.null(massir_sex_labels)){
          massir_sex_labels <- rep(NA, ncol(expr2))
        }
        massir_sex_labels[massir_sex_labels=="female"& !is.na(massir_sex_labels)] <- 0
        massir_sex_labels[massir_sex_labels=="male"& !is.na(massir_sex_labels)] <- 1
        
        #massir_lab_s  <- ifelse(massiRAcc(downsampled_df, ychr.genes)=="female", 0, 1)
        massir_lab <- c(massir_lab, massir_sex_labels)
        
        # Toker labels
        toker_keys <- keys[keys %in% toker_list2.2]
        toker_sex_labels <- tokerSexLab(expr2, f.genes=toker_list2f, m.genes=toker_list2m, min.samples=3) # /// TODO this needs to be updated
        if (is.null(toker_sex_labels)){
          toker_sex_labels <- rep(NA, ncol(expr2))
        }
        toker_sex_labels[toker_sex_labels=="female" & !is.na(toker_sex_labels)] <- 0
        toker_sex_labels[toker_sex_labels=="male" & !is.na(toker_sex_labels)] <- 1
        #toker_lab_s <- ifelse(tokerAcc(downsampled_df)=="female", 0, 1)
        
        toker_lab <- c(toker_lab, toker_sex_labels)
        
        # our sex labels
        
        downsampled_df_rank <- test_rank[,colnames(downsampled_df)]
        data(default_fit)
        rank_lab_s <- predSexLab(default_fit, downsampled_df_rank, numeric=TRUE)
        rank_lab <- c(rank_lab, rank_lab_s)
        
        actual_lab <- c(actual_lab, downsampled_labels)
      }
      
     
    }
  #}
  return(list("actual"=actual_lab, "rank"=rank_lab, "massir"=massir_lab, "toker"=toker_lab))
}

dfFracAcc <- function(frac, frac_res){
  frac_res$toker <- sapply(frac_res$toker, as.numeric)
  frac_res$massir <- sapply(frac_res$massir, as.numeric)
  
  # creates a dataframe to describe the fractional accuracy
  frac_acc_list <- list("rank"=sum(frac_res$actual==frac_res$rank, na.rm=TRUE)/length(frac_res$actual), 
                        "toker"=sum(frac_res$actual==frac_res$toker, na.rm=TRUE)/length(frac_res$actual), 
                        "toker_no_na"=sum(frac_res$actual==frac_res$toker, na.rm=TRUE)/sum(!is.na(frac_res$toker)),
                        "massir"=sum(frac_res$actual==frac_res$massir, na.rm=TRUE)/length(frac_res$actual),
                        "massir_no_na"=sum(frac_res$actual==frac_res$massir, na.rm=TRUE)/sum(!is.na(frac_res$massir)))  
  df <- data.frame(frac_acc_list)
  df <- gather(df, "method", "accuracy")
  df$frac <- frac
  return(df)
}
set.seed(1120)
#load("data/train.RData")
#load("data/extended_test.RData")
set.seed(12)
acc_res0.05 <- calcAccDown(0.05, train[[1]]$expr, train_rank, 
                          train_lab, test[[1]]$expr, test1_rank, test1_lab)
acc_res0.95 <- calcAccDown(0.95, train[[1]]$expr, train_rank, 
                           train_lab, test[[1]]$expr, test1_rank, test1_lab)


intervals <- c(0, 0.05, 0.15, 0.3, 0.6, 0.9, 0.95, 1)
acc_res <-   lapply(intervals, function(frac)
  calcAccDown(frac, train[[1]]$expr, train_rank, 
              train_lab, test[[1]]$expr, test1_rank, test1_lab))
frac_acc <- do.call(rbind, lapply(1:length(acc_res), function(i) dfFracAcc(intervals[[i]],acc_res[[i]])))
# plot the accuracy across fractions
frac_acc0.05 <- dfFracAcc(0.05, acc_res0.05)
frac_acc0.95 <- dfFracAcc(0.95, acc_res0.95)
frac_acc2 <- frac_acc %>% bind_rows(frac_acc0.05) %>% bind_rows(frac_acc0.95)
frac_acc2$method[frac_acc2$method=="rank"] <- "exprsex"
save(frac_acc2, file="data/05_performance/frac_acc2_0204.RData")
ggplot(frac_acc2, aes(frac, accuracy, col=method))+xlab("fraction males")+geom_line(size=1.5, alpha=0.9)+geom_point(size=1.5)
ggsave("figures/fractional_acc.png", dpi="print", width=5, height=3.5)
##