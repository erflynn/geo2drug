

# GOAL:
# is missingness correlated with sex label?
# are there certain genes that we should not consider missingness in?


dim(ss_rank)
ss_lab <- c(ss1_lab, ss2_lab)
ss_sm_f <- ss_rank[sl_genes$f,]
ss_sm_m <- ss_rank[intersect(sl_genes$m, rownames(ss_rank)),]
ss_stud <- read_csv("data/01_sample_lists/human_single_sex_studies.csv")


library(naniar)
#vis_miss(data.frame(cbind(ss_sm_f[,ss_lab==0], ss_sm_f[,ss_lab==1])))
#vis_miss(data.frame(cbind(ss_sm_m[,ss_lab==0], ss_sm_m[,ss_lab==1])))

count_miss <- function(my.gpl){
  gses <- ss_stud %>% filter(gpl == my.gpl) %>% select(gse) %>% unique()
  ss_lab_p <- ss_lab[filterColsByGSE(names(ss_lab), gses$gse)]
  # count missing in m vs f samples
  ss_miss_f_females <- apply(ss_sm_f[,ss_lab_p==0], 1, function(x) sum(is.na(x)))
  ss_miss_f_males <- apply(ss_sm_f[,ss_lab_p==1], 1, function(x) sum(is.na(x)))
  
  miss_f <- cbind(ss_miss_f_females, ss_miss_f_males)
  
  ss_miss_m_females <- apply(ss_sm_m[,ss_lab_p==0], 1, function(x) sum(is.na(x)))
  ss_miss_m_males <- apply(ss_sm_m[,ss_lab_p==1], 1, function(x) sum(is.na(x)))
  
  miss_m <- cbind(ss_miss_m_females, ss_miss_m_males)
  
  return(list("f"=miss_f, "m"=miss_m, 
              "num_f"=sum(ss_lab_p==0), "num_m"=sum(ss_lab_p==1)))
}
ss_gpls <- unique(ss_stud$gpl)

miss2 <- lapply(ss_gpls, count_miss)
# TODO - need to count how miss
res <- do.call(cbind, lapply(miss2, function(x) x$m))
heatmap(t(res))
res2 <- do.call(cbind, lapply(miss2, function(x) x$f))
heatmap(t(res2))
# separate this out by platform

# for each platform, is missing gene predictive of sex?