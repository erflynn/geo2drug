require('tidyverse')
source("code/utils/general_utils.R")


args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
idx <- as.numeric(args[2])
ref_dir <- sprintf("gpl_ref/%s", organism)
MIN.ROWS <- 8000

list.gpls <- read_csv(sprintf("data/01_sample_lists/%s_gpl.csv", organism))

chunk.gpls <- extractChunk(list.gpls$gpl, idx, SIZE.CHUNK=50)

lapply(chunk.gpls, function(gpl){
  print(gpl)
  tryCatch({
    ref_tab <- parse_entrez_from_gpl(gpl, ref_dir)
  }, error = function(e){
    print(sprintf("Error loading %s", gpl))
    print(e)
    ref_tab <- NULL
  })
  
  if (is.data.frame(ref_tab)){
    if (nrow(ref_tab) < MIN.ROWS){
      print(sprintf("Error, insufficient rows for %s", gpl))
    }
    save(ref_tab, file=sprintf("%s/%s_map.RData", ref_dir, gpl))
  }
})
