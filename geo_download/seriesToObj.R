
require('MetaIntegrator')
require('GEOquery')
require('RMySQL')

in_dir <- "gses/matrix"
out_dir <- "gses/rObj"
MAX_NA_FRAC <- 0.3  

options('GEOquery.inmemory.gpl'=TRUE)

checkDownload <- function(fname){
  tryCatch({
    acc <- strsplit(fname, "_")[[1]][[1]] # parse the accession from the fname
    gse <- getGEO(filename=sprintf("%s/%s", in_dir, fname)) 
    gse2 <- MetaIntegrator:::.GEOqueryGEO2GEM(gse, acc, qNorm=FALSE) # from MetaIntegrator
    
    # check comments and for lots of NAs
    if (gse2$exp_comment == "Expression data is missing"){
      cat(sprintf("%s: expression data is missing\n", acc), file=logfile, append=TRUE)
      return(NA)
    }
    key_count_nas <- table(is.na(gse2$keys))
    key_nas <- ifelse('TRUE' %in% names(key_count_nas), key_count_nas[['TRUE']]==sum(key_count_nas), FALSE )
    key_nas <- key_count_nas[['TRUE']]==sum(key_count_nas)
    if (gse2$key_comment %in% c("Annotation absent", "fData file not available") | key_nas ){
      cat(sprintf("%s: keys absent\n", acc), file=logfile, append=TRUE)
      return(NA)
    }
    exp_count_nas <- table(is.na(gse2$expr)) 
    num_cells <- sum(exp_count_nas)
    exp_na_frac <- ifelse('TRUE' %in% names(exp_count_nas), exp_count_nas[['TRUE']]/num_cells > MAX_NA_FRAC, FALSE)
    if (exp_na_frac){
      cat(sprintf("%s: over 30 percent of expression data missing\n", acc), file=logfile, append=TRUE) # <-- this is a warning
      
    }
    if (nrow(gse2$expr) < 10000){
      cat(sprintf("%s: fewer than 10k genes present on this platform\n", acc), file=logfile, append=TRUE)
    }
    #save(gse2, file=sprintf("%s/%s.RData", out_dir, acc))
    return(gse2)
  }, error = function(err){
    print(acc)
    cat(sprintf("%s unknown error during sex labeling\n", acc), file=logfile, append=TRUE)
    
    print(err)
    return(NA)
  }
  )
}


args <- commandArgs(trailingOnly=TRUE)
chunk_num <- as.numeric(args[1])

fnames <- list.files(in_dir)
logfile <- sprintf("logs/%s_%s.log", "chunk_check", chunk_num)
tmp <- chunkRun(fnames, checkDownload, out_dir, chunk_num, CHUNK_SIZE=100, GROUP_SIZE=10, log_prefix="chunk_check")



