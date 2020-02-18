
# to run:
# Rscript gse_to_mat.R gse.file out.dir idx

require('exprsex')
require('MetaIntegrator')
require('R.utils')

source("code/utils/general_utils.R")


SIZE.CHUNK <- 100


# parse arguments
args <- commandArgs(trailingOnly=TRUE)
gse.file <- args[1]
gse.list <- read.csv(gse.file, header=TRUE, stringsAsFactors=FALSE)
OUT.DIR <- args[2]

organism <- args[3]
GSE.DIR <- sprintf("gses/%s/matrix/", organism)
GPL.DIR <- sprintf("gpl_ref/%s/", organism)

idx <- as.numeric(args[4])
print(idx)


list.gses <- gse.list[,1]

chunk.gses <- extractChunk(list.gses, idx, SIZE.CHUNK)
print(chunk.gses)


downloadData <- function(gse){
  print(gse)
  
  out.file <- sprintf("%s/%s_mat.RData", OUT.DIR, gse)
  if (!file.exists(out.file)){
    geo.obj <- NULL
    # try to download the file
    try(geo.obj <- exprsex::getPrepGSE(gse, gpl.dir=GPL.DIR, gse.dir=GSE.DIR)) 
    if (is.null(geo.obj)){
      print(sprintf("Could not load %s", gse))
    } else {
      save(geo.obj, file=out.file)
      print(sprintf("Succeeded with %s", gse))
      
    }
  } else {
    print(sprintf("Already loaded %s", gse))
  }
}


gses.to.run <- setdiff(chunk.gses, c("GSE25219", "GSE18927", "GSE28387", "GSE30727", "GSE37138", "GSE50421", "GSE77714", "GSE19090"))

lapply(gses.to.run, function(gse){
  tryCatch({
  res <- withTimeout({
    downloadData(gse)
    }, timeout=18000)
  }, TimeoutException = function(ex) {
    message(sprintf("Timeout. Skipping %s.", gse))
  })
})