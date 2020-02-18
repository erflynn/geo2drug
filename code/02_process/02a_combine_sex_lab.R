
# combine the sex labels together

require('tidyverse')
source("code/utils/general_utils.R")

SIZE.CHUNK <- 100
args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
ds.dir <- args[2] # the dataset directory (e.g. '05_single_sex' or '03_silver_std')
filter.tab <- as.logical(args[3]) # whether to filter the file based on annot
idx <- as.numeric(args[4])

out.dir <- sprintf("data/%s/%s/02_keep_labels/", ds.dir, organism)

all_f <- list.files(path=sprintf("data/%s/%s/01_compare_labels/", ds.dir, organism))
gse.list <- extractChunk(all_f, idx, SIZE.CHUNK)


lapply(gse.list, function(f){
  load(sprintf("data/%s/%s/01_compare_labels/%s",ds.dir, organism,f)) # --> pheno2
  gse.id <-strsplit(f, "_")[[1]][[1]]
  pheno_tab <- pheno2[,c("gsm", "toker_sex", "massir_sex", "text_sex")]
  if (!filter.tab){
    write_csv(pheno_tab, sprintf("%s/%s_pheno.csv", out.dir, gse.id))
  } else {

    if (organism == "human" & !all(is.na(pheno_tab$toker_sex))){
      matched_lab <- pheno_tab %>% dplyr::filter(toker_sex==massir_sex & massir_sex==text_sex)	
    } else {
      print(table(pheno_tab[,c("massir_sex", "text_sex")]))
      matched_lab <- pheno_tab %>% dplyr::filter(massir_sex==text_sex)
    }
    counts <- matched_lab %>% filter(text_sex %in% c("male", "female")) %>% group_by(text_sex) %>% count()
    print(counts)
    if(nrow(matched_lab) > 1 & all(counts$n >= 10)){ 
      # // TODO: BUGGGGGG!!! if none of something it shouldn't be here
    	print(gse.id)
    	write_csv(matched_lab, sprintf("%s/%s_pheno.csv", out.dir, gse.id))
    } else {
    	 print(sprintf("%s nsufficient after filtering", gse.id))
    }
  }
})

