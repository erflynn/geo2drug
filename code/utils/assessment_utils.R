# ---- helper functions ---- #
read_process_mat <- function(path){
  df <- data.table::fread(path, data.table=FALSE)
  rownames(df) <- df$V1
  df$V1 <- NULL
  return(as.matrix(df))
}

read_process_lab <- function(path){
  lab <- read.csv(path)
  lab <- lab %>% 
    mutate("sample"=paste(gse, gsm, sep=".")) %>% 
    mutate(sample=stringr::str_replace_all(sample, "-", "\\."))
  sex_lab <- lab$sex
  names(sex_lab) <- lab$sample
  return(sex_lab)
}

read_process_ds <- function(organism, sub.dir, prefix){
  dir.path <- sprintf("data/%s/%s/03_out_mat/",sub.dir, organism)
  my.f <- list.files(path=dir.path, pattern=sprintf("dat_%s_expr", prefix))
  suffix.list <- unlist(
    stringr::str_extract_all(stringr::str_extract_all( my.f, "expr_[0-9]+.csv"),"[0-9]+"))
  # pre-process
  loaded.dat <- lapply(suffix.list, function(idx){
    dat_expr <- read_process_mat(
      sprintf("%s/dat_%s_expr_%s.csv", dir.path, prefix, idx))
    dat_lab <- read_process_lab(
      sprintf("%s/dat_%s_lab_%s.csv", dir.path, prefix, idx))
    overlap.ids <- intersect(names(dat_lab), colnames(dat_expr))
    
    return(list("expr"=dat_expr[,overlap.ids], "lab"=dat_lab[overlap.ids]))
    
  })
  return(loaded.dat)
}

# // TODO: deprecated
# readInTest <- function(organism, idx, test_type, genes){
#   # read in the rest of the test data #
#   test_rank <- read_process_mat(
#     sprintf("data/03_silver_std/%s/03_out_mat/test_%s_rank_%s.csv", organism, test_type, idx), genes)
#   test_lab <- read_process_lab(sprintf("data/03_silver_std/%s/03_out_mat/test_%s_lab_%s.csv", organism, test_type, idx))
#   
#   return(list("rank"=test_rank[,names(test_lab)], "lab"=test_lab))
# }

read_clean_sl_genes <- function(organism){
  fgenes_df <- read.csv(sprintf("data/03_silver_std/%s/04_meta_res/fgenes.csv", organism))
  mgenes_df <- read.csv(sprintf("data/03_silver_std/%s/04_meta_res/mgenes.csv", organism))
  fgenes_df2 <- fgenes_df %>% 
    mutate(id=as.character(entrezgene_id)) %>% 
    select(id) %>%
    unique() 
  mgenes_df2 <- mgenes_df %>% 
    mutate(id=as.character(entrezgene_id)) %>% 
    select(id) %>%
    unique() 
  return(list("f"=fgenes_df2$id, "m"=mgenes_df2$id))
}

confMat <- function(actual, pred){
  table(data.frame(cbind(actual, pred)))
}
calcAcc <- function(actual, pred){
  sum(pred==actual)/length(actual)
}


getGSEsDotted <- function(list.names){
  sapply(list.names, function(x) {
    z <- strsplit(x, "\\.")[[1]]; ifelse(length(z)==2, z[[1]], paste(z[[1]], z[[2]], sep="-"))})
}

sepReformatGPL <- function(df) {
  df_s <- df %>% filter(!str_detect(gpl, "\\|")) 
  df_m <- df %>% filter(str_detect(gpl, "\\|"))  %>% separate_rows(gpl, sep="\\|") %>%
    mutate(gse=sprintf("%s-%s", gse, gpl)) 
  return(data.frame(rbind(df_s, df_m)))
}


pred_out_df <- function(fit, df, lab, organism, ref){
  pred_df <- predSexLab(fit, df, scores=TRUE)
  pred_df <- pred_df %>% rename(pred_sex=sex) 
  pred_df$true_sex <- lab
  pred_df$gsm <- sapply(rownames(pred_df ), function(x) {l <- strsplit(x, "\\.")[[1]]; l[[length(l)]]})
  pred_df$gse <- getGSEsDotted(rownames(pred_df))
  
  pred_gpl <- left_join(pred_df, ref %>% select(gse, gpl) %>% unique())
  return(pred_gpl)
}


sampWrapper <- function(df, NUM_SEL=4){
  gpl_sm <- df %>% group_by(gpl) %>% summarize(n=length(gse)) %>% filter(n < NUM_SEL)
  df_down <- df %>% filter(! (gpl %in% gpl_sm$gpl)) %>% group_by(gpl) %>% sample_n(NUM_SEL) %>% ungroup()
  df_sm <- df %>% filter(gpl %in% gpl_sm$gpl)
  return(data.frame(rbind(df_down, df_sm)))
}

acc_by_study <- function(pred_df){
  pred_df %>% 
    group_by(gpl, gse) %>%
    summarize(true_pred=sum(true_sex==pred_sex, na.rm=TRUE), 
              num_samples=length(true_sex)) %>%
    mutate(acc=true_pred/num_samples) %>%
    ungroup()
}


summarizeAcc <- function(pred_df, lab, organism){
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
  write_csv(pred_acc2, sprintf("data/05_performance/%s_pred_acc_by_gpl_%s.csv", organism, lab))
  
  ## we see that much of this *IS* separable using the score metric, 
  ##    but there is large platform-to-platform variability
  ggplot(pred_df, aes(x=score_f, y=score_m))+geom_point(aes(color=true_sex, alpha=0.3))+
    geom_abline(slope=1, intercept=fit$threshold)+facet_wrap(vars(gpl))
  ggsave(sprintf("figures/plat_performance/%s_by_plat_%s.png", organism, lab), 
         dpi="print", width=11, height=10)
  
  # overall
  ggplot(pred_df, aes(x=score_f, y=score_m))+geom_point(aes(color=true_sex, alpha=0.3))+
    geom_abline(slope=1, intercept=fit$threshold) # TODO - add platform here
  ggsave(sprintf("figures/plat_performance/%s_overall_%s.png", organism, lab), 
         dpi="print", width=11, height=10)
}
