
# combine the sex labels together

require('tidyverse')




SIZE.CHUNK <- 100
args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
idx <- as.numeric(args[2])

out.dir <- sprintf("data/03_silver_std/%s/02_keep_labels/", organism)

all_f <- list.files(path=sprintf("data/03_silver_std/%s/01_compare_labels/", organism))
NUM.CHUNKS <- ceiling(length(all_f)/SIZE.CHUNK)
end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(all_f), (idx+1)*SIZE.CHUNK)

gse.list <- all_f[(idx*SIZE.CHUNK+1):end_idx]


lapply(gse.list, function(f){
  load(sprintf("data/03_silver_std/%s/01_compare_labels/%s",organism,f)) # --> pheno2
  gse.id <-strsplit(f, "_")[[1]][[1]]
  pheno_tab <- pheno2[,c("gsm", "toker_sex", "massir_sex", "text_sex")]
  if (organism == "human"){
    matched_lab <- pheno_tab %>% dplyr::filter(toker_sex==massir_sex & massir_sex==text_sex)
  } else {
    print(table(pheno_tab[,c("massir_sex", "text_sex")]))
    matched_lab <- pheno_tab %>% dplyr::filter(massir_sex==text_sex)
  }
  counts <- matched_lab %>% filter(text_sex %in% c("male", "female")) %>% group_by(text_sex) %>% count()
  print(counts)
  if(nrow(matched_lab) > 1 & all(counts$n >= 10)){
    print(gse.id)
    write_csv(matched_lab, sprintf("%s/%s_pheno.csv", out.dir, gse.id))
  } else {
    print(sprintf("%s nsufficient after filtering", gse.id))
  }
})

