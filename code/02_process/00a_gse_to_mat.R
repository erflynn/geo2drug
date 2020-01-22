
# to run:
# Rscript gse_to_mat.R gse.file out.dir idx

require('exprsex')
require('MetaIntegrator')
require('R.utils')
SIZE.CHUNK <- 50

# parse arguments
args <- commandArgs(trailingOnly=TRUE)
gse.file <- args[1]
gse.list <- read.csv(gse.file, header=TRUE)
OUT.DIR <- args[2]

organism <- args[3]
GSE.DIR <- sprintf("gses/%s/matrix/", organism)
GPL.DIR <- sprintf("gpl_ref/%s/", organism)

idx <- as.numeric(args[4])
print(idx)



NUM.CHUNKS <- ceiling(nrow(gse.list)/SIZE.CHUNK)
end_idx <- ifelse((NUM.CHUNKS-1) == idx ,nrow(gse.list), (idx+1)*SIZE.CHUNK)
gse.list <- gse.list[(idx*SIZE.CHUNK+1):end_idx,]

print(gse.list[,1])

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

lapply(gse.list[,1], function(gse){
  tryCatch({
  res <- withTimeout({
    downloadData(gse)
    }, timeout=18000)
  }, TimeoutException = function(ex) {
    message(sprintf("Timeout. Skipping %s.", gse))
  })
})