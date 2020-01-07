# create_exp_mat.R
# E Flynn
#
# Code for creating an expression matrix from the data

require('miceadds')
require('MetaIntegrator')
require('tidyverse')

PLATFORM_DIR <- "gpls/gpl_mapping"

#' Grab the gene to probe mapping for a given object
#'
#' @param gse.obj - the gse object of interest
#' @return gene.to.probe - a mapping from genes to probes
getGeneToProbe <- function(gse.obj){
  
  # check if the mapping already exists - load if it does
  platform <- gse.obj$platform
  platform.path <- sprintf("%s/%s.RData", PLATFORM_DIR, platform)
  if (file.exists(platform.path)){
    load.Rdata(platform.path, "gene.to.probe")
  }  else {
    # extract the keys and then convert it to a table
    keys <- gse.obj$keys
    list.keys <- keys
    key.df <- data.frame(list.keys, names(list.keys))
    colnames(key.df) <- c("gene", "probes")
    key.df$gene <- as.character(key.df$gene)
    key.df$probes <- as.character(key.df$probes)
    key.df2 <- separate_rows(key.df, gene, sep=",") # split cases where probe --> mult genes
    gene.to.probe <- split(key.df2$probes,  key.df2$gene)
    save(gene.to.probe, file=platform.path)
  }
  return(gene.to.probe)
}


#' Summarizes MetaIntegrator to the gene level
#'
#' @param gse.obj
#' @param gene.list
#'
#' @return expression matrix
convertToGenes <- function(gse.obj, gene.list){
  expData <- gse.obj$expr
  gene.to.probe <- getGeneToProbe(gse.obj)
  gene.to.probe <- gene.to.probe[(names(gene.to.probe) %in% gene.list)]
  expData2 <- do.call(cbind, lapply(1:length(gene.to.probe), function(x) {
    # get the gene and the probe
    g <- names(gene.to.probe)[x]
    p <- unlist(gene.to.probe[g])
    if (length(p)>1){ 
      expD <- expData[p,]
      df <- (data.frame(colMeans(expD, na.rm=TRUE)))
      return(df)
    }
    else {
      df <- data.frame(expData[p,])
      return(df)
    }})) ### ALSO SLOW...
  
  colnames(expData2) <- names(gene.to.probe)
  expData2.2 <- data.frame(t(expData2)) # columns are samples, rows are genes
  
  # create a data fram of NAs for missing genes    
  missing.genes <- setdiff(gene.list, list.keys)
  missing.vec <- rep(NA, ncol(expData2.2))
  missing.df <- do.call(rbind, lapply(1:length(missing.genes), function(x) missing.vec))
  rownames(missing.df) <- missing.genes
  colnames(missing.df) <- colnames(expData2.2)
  
  # put together and reorder
  expDataPlusMiss <- rbind(expData2.2, missing.df )
  expData2.3 <- expDataPlusMiss[gene.list,] # REORDER so it matches other data
  
  return(expData2.3)
}

# summarize each dataset to an expression matrix
combineToExpMat <- function(obj.files){
  # convert each to genes
  expMat2 <- do.call(cbind,
                     lapply(obj.files, function(obj.file){
                       load.Rdata( sprintf("%s/%s", in_dir, obj.file), "chunk_ds")
                       
                       gse_obj_list <- chunk_ds[!is.na(chunk_ds)]   # skip the ones that are NAs
                       expMat <- do.call(cbind, 
                                         lapply(gse_obj_list, function(gse_obj) convertToGenes(gse_obj, list_genes))) ## TODO: parallelize
                       return(expMat)
                     }))
  #return(expMat2)
  write.csv(expMat2, file=sprintf("%s/expData_%s.csv", out_dir, chunk_num), row.names=TRUE, quote=FALSE)
}

args <- commandArgs(trailingOnly=TRUE)
chunk_num <- as.numeric(args[1])
in_dir <- spritnf("gses_%s/rObj", dir_id)

list_genes <- read.csv("GPL570_genes.csv")[,1]
out_dir <- sprintf("gses_%s/exp_mat", dir_id)

list.items <- list.files(in_dir, sprintf("chunk%s_", chunk_num))
combineToExpMat(list.items)

#tmp <- chunkRun(obj.files, combineToExpMat, out_dir, chunk_num, CHUNK_SIZE=100, GROUP_SIZE=10, log_prefix="chunk_expr")
