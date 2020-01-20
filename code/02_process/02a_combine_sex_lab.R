


require('tidyverse')
out.dir <- "/scratch/users/erflynn/sex_labeling/geo_pipeline/data/keep_labels/"



SIZE.CHUNK <- 100
args <- commandArgs(trailingOnly=TRUE)
idx <- as.numeric(args[1])

all_f <- list.files()
NUM.CHUNKS <- ceiling(length(all_f)/SIZE.CHUNK)
end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(all_f), (idx+1)*SIZE.CHUNK)

gse.list <- all_f[(idx*SIZE.CHUNK):end_idx]


lapply(gse.list, function(f){
  load(f) # --> pheno2
  gse.id <-strsplit(f, "_")[[1]][[1]]
  pheno_tab <- pheno2[,c("gsm", "toker_sex", "massir_sex", "text_sex")]
  matched_lab <- pheno_tab %>% dplyr::filter(toker_sex==massir_sex & massir_sex==text_sex)
  counts <- matched_lab %>% filter(text_sex %in% c("male", "female")) %>% group_by(text_sex) %>% count()
  print(counts)
  if(nrow(matched_lab) > 1 & all(counts$n >= 10)){
    print(gse.id)
    write_csv(matched_lab, sprintf("%s/%s_pheno.csv", out.dir, gse.id))
  } else {
    print(sprintf("%s nsufficient after filtering", gse.id))
  }
})

