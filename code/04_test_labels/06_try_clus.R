# Goal: see if the comparison of exprsex + duda-hart test improves performance
#  this is of particular concern for "rescuing" single-sex studies

# run a duda-hardt test


# 1) duda-hart test --> does it work?
require('fpc')
require('cluster')

# how does it perform on single-sex studies?


run_all <- function(gse, my.df){
  keep.cols <- filterColsByGSE(colnames(my.df), c(gse))
  gse1 <- my.df[,keep.cols]
  sl_genes <- c(sex_lab_genes$f, sex_lab_genes$m)
  expMat <- gse1[sl_genes,]
  expMat2 <- expMat[complete.cases(expMat),]
  if(nrow(expMat2) < 2 | ncol(expMat) < 2){
    dh <- NA
    clus.widths <- NA
    tL <- NA
  } else {
    #clusT <- kmeans(t(expMat2), 2)
    #table(clusT$cluster)
    pr <- pam(t(expMat2), k=2)
    clus.widths <- pr$silinfo$clus.avg.widths # average widths  
    dh <- dudahart2(t(expMat2), pr$clustering, alpha=0.001)
    # cluster label
    # filter -OR- impute
    tL <- tokerSexLab(expMat2, f.genes=sex_lab_genes$f, m.genes=sex_lab_genes$m)
    # look at how it corresponds to exprsex label
  }
  
  eL <- predSexLab(fit, gse1, numeric_lab=FALSE)
  df <- list("cw"=clus.widths, "dh"=dh, "el"=eL, "tL"=tL)
  return(df)
}

test_duda <- lapply(unique(test_pred_df2$gse), function(gse){
  print(gse)
  return(run_all(gse, test_rank))
})


ss_duda <- lapply(unique(ss_pred_df$gse), function(gse){
  print(gse)
  if (gse %in% c("GSE45898", "GSE66667", "GSE68021")){
    return(NA)
  }
  return(run_all(gse, ss_rank))
})

ss_duda_df <- data.frame(do.call(rbind, ss_duda[!is.na(ss_duda)]))

lapply(test_duda, function(x) ifelse(length(x$cw)==1, NA, x$dh$cluster1))
t(sapply(test_duda, function(x) x$cw))

table(sapply(test_duda, function(x) any(x$tL!=x$el))) 
# 39 of these have a mislabeled one
table() 
# 22 have a disagreement
table(sapply(ss_duda[!is.na(ss_duda)], function(x) x$dh$cluster1))
# 18 are TRUE
comb_ss <- 
  do.call(rbind,
          list(
               sapply(ss_duda[!is.na(ss_duda)], function(x) x$dh$cluster1),
        sapply(ss_duda[!is.na(ss_duda)], function(x) x$cw),
        sapply(ss_duda[!is.na(ss_duda)], function(x) length(x$tL)==0),
        sapply(ss_duda[!is.na(ss_duda)], function(x) 1-sum(x$tL!=x$el)/length(x$el))
        ))

colnames(comb_ss) <- unique(ss_pred_df$gse)[!is.na(ss_duda)]
comb_df <- data.frame(t(comb_ss))
colnames(comb_df) <- c("dh", "cw1", "cw2", "no_toker", "acc")
comb_df %>% filter(dh==1)
comb_df %>% filter(no_toker==0 & dh==0) %>% arrange(desc(acc))

test_duda_df <- data.frame(do.call(rbind, test_duda[!is.na(test_duda)]))

# ranked in increasing order --> 
#  lower rank == lower value
# so HIGHER value of genes --> HIGHER rank
#  higher M than f, more likely M than F
keep.cols <- filterColsByGSE(colnames(ss[[1]]$expr), c("GSE41790"))
gse1 <- ss[[1]]$expr[,keep.cols]
sl_genes <- c(sex_lab_genes$f, sex_lab_genes$m)
expMat <- gse1[sl_genes,]
median(apply(expMat[sex_lab_genes$m, ], 2, median, na.rm=TRUE))
median(apply(expMat[sex_lab_genes$f, ], 2, median, na.rm=TRUE))


# PCA - x
miceadds::load.Rdata(sprintf("../gpl_ref/%s_xy_genes.RData", organism), "xy_genes") 
xy_genes2 <- xy_genes %>% unique() %>% filter(!is.na(gene))
x_genes <- xy_genes2 %>% filter(chr=="X") %>% dplyr::select(gene)
y_genes <- xy_genes2 %>% filter(chr=="Y") %>% dplyr::select(gene)
y2 <- intersect(y_genes$gene, rownames(test[[1]]$expr))
x2 <- intersect(x_genes$gene, rownames(test[[1]]$expr))
comb.df <-cbind(test[[1]]$expr, test[[2]]$expr)
comb.lab <-c(test[[1]]$lab, test[[2]]$lab)
comb.df[is.na(comb.df)] <- 0
### ok so looks like we need some imputation....

x3 <- x2[apply(comb.df[x2,], 1, function(x) sum(is.na(x)))==0]

y3 <- y2[apply(comb.df[y2,], 1, function(x) sum(is.na(x)))==0]

pcy <- prcomp(comb.df[sex_lab_genes$m,])
pcx <- prcomp(comb.df[sex_lab_genes$f,])
plot(pcx$rotation[,1], 
     pcy$rotation[,1], 
     col=ifelse(comb.lab==0, "red", "blue"))

