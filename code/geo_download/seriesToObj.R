# seriesToObj.R
# E Flynn
#
# Code for converting a series matrix to a MetaIntegrator R object.
# In the process, it checks for missing expression data, >30% NAs

require('MetaIntegrator')
require('GEOquery')
require('RMySQL')


MAX_NA_FRAC <- 0.3  
MIN_GENES <- 10000

options('GEOquery.inmemory.gpl'=TRUE)

checkDownload <- function(fname){
  tryCatch({
    acc <- strsplit(fname, "_")[[1]][[1]] # parse the accession from the fname
    print(acc)
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
      cat(sprintf("%s: over %d percent of expression data missing\n", acc, MAX_NA_FRAC*100), file=logfile, append=TRUE) # <-- this is a warning
    }
    #save(gse2, file=sprintf("%s/%s.RData", out_dir, acc))
    if (nrow(gse2$expr) < MIN_GENES){
      cat(sprintf("%s: fewer than %s genes present on this platform\n", acc, MIN_GENES), file=logfile, append=TRUE)
    }
    return(gse2)
  }, error = function(err){
    print(acc)
    cat(sprintf("%s unknown error during conversion to obj\n", acc), file=logfile, append=TRUE)
    
    print(err)
    return(NA)
  }
  )
}


args <- commandArgs(trailingOnly=TRUE)
chunk_num <- as.numeric(args[1])
dir_id <- args[2]

in_dir <- sprintf("gses_%s/matrix", dir_id)
#out_dir <- sprintf("gses_%s/rObj", dir_id)


fnames <- list.files(in_dir)
list.gses <- fnames
logfile <- sprintf("logs_%s/check_chunk%s.log", dir_id, chunk_num)
cat(" \n", file=logfile)

CHUNK_SIZE=100
if (chunk_num*CHUNK_SIZE > length(list.gses)){
  chunk_gses <- list.gses[((chunk_num-1)*CHUNK_SIZE):length(list.gses)]
} else {
  chunk_gses <- list.gses[((chunk_num-1)*CHUNK_SIZE):(chunk_num*CHUNK_SIZE)]
}
print(chunk_gses)

GROUP_SIZE = 10
NUM_GROUPS= ceiling((length(chunk_gses))/GROUP_SIZE)

# run in chunks of ten, then save
for (i in 1:NUM_GROUPS){
  tryCatch({
  if (i==NUM_GROUPS){
    chunk_gses_i <- chunk_gses[((i-1)*GROUP_SIZE+1):length(chunk_gses)]
  } else {
    chunk_gses_i <- chunk_gses[((i-1)*GROUP_SIZE+1):(i*GROUP_SIZE)]
  }
  chunk_gses_i <- chunk_gses_i[!is.na(chunk_gses_i)]
  chunk_ds <- lapply(chunk_gses_i, checkDownload)
  names(chunk_ds) <- chunk_gses_i
  save(chunk_ds, file=sprintf("gses_%s/rObj/chunk%s_%s.RData", dir_id, chunk_num, i))
}, error = function(err){
    print(err)
})
}

