
# code for redownloading GPLs that failed


library(parallel)
library(tictoc)

human_gpls <- read_csv("gpls/human_gpl.csv")
require('profvis')


list_downloaded <- read_table("gpls/list_human_downloaded.txt", skip=1, col_names=FALSE)
gpl_downloaded <- sapply(list_downloaded$X9, function(x) strsplit(x, ".", fixed=TRUE)[[1]][[1]])
gpl_to_download <- setdiff(human_gpls$gpl, gpl_downloaded)

downloadMapGPL <- function(gpl){
  tryCatch({
    tmp <- exprsex:::.getGeneToProbe(gpl, platform_dir="gpls/gpl_download/");
    return(sprintf("%s succeeded", gpl))
  }, error=function(err) {
    return(sprintf("%s failed: %s", gpl,err))
  })
}


profvis({
  cl <- makeCluster(detectCores() - 1)
  output <- 
    parSapply( cl, gpl_to_download[51:200], downloadMapGPL)
  stopCluster(cl)
})
output1 <- output 
failures1 <- names(output1)[sapply(output1, str_length) > 20]



# -- this fails... too much memory?? -- #



library(rslurm)
require('exprsex')
require('MetaIntegrator')
require('tidyverse')
downloadMapGPL <- function(gpl){
  tryCatch({
    tmp <- exprsex:::.getGeneToProbe(gpl, platform_dir="/scratch/users/erflynn/sex_labeling/download_geo_v2/gpls/gpl_download/");
    return(sprintf("%s succeeded", gpl))
  }, error=function(err) {
    return(sprintf("%s failed: %s", gpl,err))
  })
}

human_gpls <- read_csv("gpls/human_gpl.csv")
list_downloaded <- read_table("gpls/list_human_downloaded.txt", skip=1, col_names=FALSE)
gpl_downloaded <- sapply(list_downloaded$X9, function(x) strsplit(x, ".", fixed=TRUE)[[1]][[1]])
gpl_to_download <- setdiff(human_gpls$gpl, gpl_downloaded)

sopt <- list(partition = 'rbaltman', error="gpl%A_%a.err")
pars <- data.frame("gpl"=gpl_to_download)

sjob <- slurm_apply(downloadMapGPL, pars, jobname = 'download_annot',
                    nodes = 60, cpus_per_node = 1, submit = TRUE, 
                    slurm_options=sopt)
