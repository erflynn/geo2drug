
require('miceadds')

# summarize each dataset to an expression matrix
combineToExpMat <- function(obj.files){
  # convert each to genes
  expMat2 <- do.call(cbind,
                     lapply(obj.files, function(obj.file){
                       load.Rdata( sprintf("%s/%s", in_dir, ßfname), "chunk_ds")
                       
                       load(obj.file)
                       gse_obj_list <- chunk_ds[!is.na(gse_obj_list)]   # skip the ones that are NAs
                       expMat <- do.call(cbind, 
                                         lapply(gse_obj_list, function(gse_obj) MetaIntegrator::getSampleLevelGeneData(gse_obj, list_genes))) ## TODO: parallelize
                       return(expMat)
                     }))
  #return(expMat2)
  write.csv(expMat2, file=sprintf("%s/expData_%s.csv", out_dir, chunk_num), row.names=, quote=FALSE)
}

args <- commandArgs(trailingOnly=TRUE)
chunk_num <- as.numeric(args[1])
list.items <- list.files(in_dir)
in_dir <- "gses/check_sex_lab"
obj.files <- list.files(in_dir)
list_genes <- "" ### READ THIS IN
out_dir <- "gses/expr_mat"

tmp <- chunkRun(obj.files, combineToExpMat, out_dir, chunk_num, CHUNK_SIZE=100, GROUP_SIZE=10, log_prefix="chunk_expr")
