# study_sex_label.R
# E Flynn
# 
# Code for running other sex labeling methods

require('tidyverse')
require('MetaIntegrator')
require('massiR')
require('miceadds')

source("code/utils/sex_lab_utils.R")



options(stringsAsFactors=FALSE)

SIZE.CHUNK <- 50

# parse arguments
args <- commandArgs(trailingOnly=TRUE)

OUT.DIR <- args[1]
MAT.DIR <- args[2]
organism <- args[3]
idx <- as.numeric(args[4])
print(idx)
logfile <- "tmp.log"
out_dir <- OUT.DIR


gse.list <- sapply(list.files(MAT.DIR), function(x) strsplit(x, "_")[[1]][[1]])

NUM.CHUNKS <- ceiling(length(gse.list)/SIZE.CHUNK)
end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(gse.list), (idx+1)*SIZE.CHUNK)
gse.list <- gse.list[(idx*SIZE.CHUNK+1):end_idx]

print(gse.list)

combined_data2 <- read.csv(sprintf("data/01_sample_lists/%s_ale_sex_lab.csv", organism))
larger_annot2 <- combined_data2
larger_annot2$text_sex <- sapply(larger_annot2$Gender, function(x)
  ifelse(x=="M", "male", ifelse(x=="F", "female", x)))

miceadds::load.Rdata(sprintf("gpl_ref/%s_xy_genes.RData", organism), "xy_genes") 
xy_genes2 <- xy_genes %>% unique() %>% filter(!is.na(gene))
x_genes <- xy_genes2 %>% filter(chr=="X") %>% dplyr::select(gene)
y_genes <- xy_genes2 %>% filter(chr=="Y") %>% dplyr::select(gene)

xy_genes2.2 <- xy_genes2$gene

# we only do toker labeling for human 
if (organism == "human"){
   miceadds::load.Rdata("gpl_ref/human_gene_map.RData", "gene_map")
   toker_list <- c("XIST", "KDM5D", "RPS4Y1")
   toker_list2 <- gene_map %>% filter(hgnc_symbol %in% toker_list) %>% dplyr::select(hgnc_symbol, entrezgene_id) %>% unique()
   toker_list2.2 <- toker_list2$entrezgene_id
   toker_list2f <- (toker_list2 %>% filter(hgnc_symbol=="XIST") )$entrezgene_id
   toker_list2m <- (toker_list2 %>% filter(hgnc_symbol %in% c("KDM5D", "RPS4Y1") ))$entrezgene_id
}

NUM_COUNTS <- 10
NA.MAX <- 0.3



studySexLabel <- function(geo.obj, gse.id){

  # // TODO - need to iterate thru
  res <- tryCatch({
    # for a GSE object, sex label and figure out if the sex labels match
    # find the ones that actually match the GSE
    to.keep <- stringr::str_detect(names(geo.obj), sprintf("%s[-|_]+", gse.id))
    geo.obj2 <- geo.obj[to.keep]
    if (length(geo.obj2)==2){
      print(sprintf("missing data for %s", gse.id))
      return(NA)
    }
    res2 <- lapply(1:length(geo.obj2), function(i){
    	 gse.obj <- geo.obj2[[i]]
      pheno <- gse.obj$pheno
      
      # filter the gene matrix to remove genes with more than 30% NAs
      expr <- gse.obj$gene_mat
      
      na.counts <- apply(expr, 1, function(x) sum(is.na(x)))
      expr2 <- expr[floor(na.counts/ncol(expr)) <= NA.MAX,]
      keys <- rownames(expr2)
      list.keys <- keys[keys %in% xy_genes2.2]
      
      if (length(list.keys) < 1){
        print(sprintf("%s missing keys", gse.id))
        cat(sprintf("%s missing keys\n", gse.id), file=logfile, append=TRUE)
        return(NA)
      }
      print("Succeeded at download")	
      ychr_keys <- keys[keys %in% y_genes$gene]
      if (organism == "human"){
      	 toker_keys <- keys[keys %in% toker_list2.2]
	       toker_sex_labels <- tokerSexLab(expr2, f.genes=toker_list2f, m.genes=toker_list2m) # /// TODO this needs to be updated
      } else {
      	toker_sex_labels <- c()
      }
      
      
      massir_sex_labels <- massiRAcc(expr2, ychr_keys) # // this works!
      
      if(length(toker_sex_labels)==0){
        cat(sprintf("%s Toker sex labeling failed\n", gse.id), file=logfile, append=TRUE)
        pheno$toker_sex <- NA
        toker_failed <- TRUE
      } else {
        toker_failed <- FALSE
        pheno$toker_sex <- toker_sex_labels
      }
      if(length(massir_sex_labels)==0){ # UPDATE
        cat(sprintf("%s MassiR sex labeling failed\n", gse.id), file=logfile, append=TRUE)
        massir_failed <- TRUE
        pheno$massir_sex <- NA
      } else {
        massir_failed <- FALSE
        pheno$massir_sex <- massir_sex_labels
      }
      if (toker_failed & massir_failed){
        cat(sprintf("%s both expression-based sex labeling failed, please skip this dataset\n", gse.id), file=logfile, append=TRUE)
        
        print(sprintf("%s Error - both expression-based sex labeling failed. please skip this dataset.", gse.id))
        return(NA)
        
      } else {
        pheno$gsm <- rownames(pheno)
        pheno2 <- left_join(pheno, larger_annot2[,c("gsm", "text_sex")], by="gsm")
        print(table(pheno2[,c("text_sex", "toker_sex", "massir_sex")]))
	save(pheno2, file=sprintf("%s/%s_sex_lab.RData", OUT.DIR,names(geo.obj)[[i]]))
	return(1)
      }
    })
  }, error = function(err){ # catch loading errors
    print(sprintf("%s error loading", gse.id))
    cat(sprintf("%s unknown error during sex labeling\n", gse.id), file=logfile, append=TRUE)
    print(err)
    return(NA)
  })
}




lapply(gse.list, function(gse.id){
print(gse.id)
gse.f <- sprintf("%s/%s_mat.RData", MAT.DIR, gse.id)

if (file.exists(gse.f)){
		     miceadds::load.Rdata(gse.f, "geo.obj")
		     res <- studySexLabel(geo.obj, gse.id)
}
		 
})