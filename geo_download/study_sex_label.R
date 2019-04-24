require('tidyverse')
require('MetaIntegrator')
require('massiR')
require('miceadds')

source("../runDownloadLabel/sex_lab_utils.R")

options(stringsAsFactors=FALSE)


combined_data2 <- read.csv("../runDownloadLabel/data/ale_combined_data.csv") ##
larger_annot2 <- combined_data2
larger_annot2$text_sex <- sapply(larger_annot2$text_sex, function(x)
  ifelse(x=="M", "male", ifelse(x=="F", "female", x)))

y_genes <- read.csv("../runDownloadLabel/data/ychr_genes.csv")
x_genes <- read.csv("../runDownloadLabel/data/xchr_genes.csv")
xy_genes <- c(x_genes[,"x"], y_genes[,"x"])
toker_list <- c("XIST", "KDM5D", "RPS4Y1")
NUM_COUNTS = 10


studySexLabel <- function(gse.obj){
  # for a GSE object, sex label and figure out if the sex labels match
  gse.id <- names(gse.obj)
  res <- tryCatch({
    gse.id <- gse.obj$name
    
    pheno <- gse.obj$pheno
    keys <- gse.obj$keys
    list.keys <- keys[keys %in% xy_genes]
    
    if (length(list.keys) < 1){
      print(sprintf("%s missing keys", gse.id))
      cat(sprintf("%s missing keys\n", gse.id), file=logfile, append=TRUE)
      return(NA)
    }
    print("Succeeded at download")	
    ychr_keys <- keys[keys %in% y_genes[,"x"]]
    toker_keys <- keys[keys %in% toker_list]
    toker_sex_labels <- imputeSex(gse.obj)
    exp <- gse.obj$expr
    massir_sex_labels <- massiRAcc(exp, names(ychr_keys))
    
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
      pheno2 <- left_join(pheno, larger_annot2[,c("gsm", "text_sex", "expr_sex")], by="gsm")
      print(table(pheno2[,c("text_sex", "toker_sex", "massir_sex")]))
      
      # if two methods --> use ones where two agree
      if (toker_failed){
        filt_pheno <- pheno2[pheno2$text_sex == pheno2$massir_sex,]
      } else {
        if (massir_failed){
          filt_pheno <- pheno2[pheno2$text_sex == pheno2$toker_sex,]
        } else {       # if all three methods --> use ones where all agree
          filt_pheno <- pheno2[(pheno2$text_sex == pheno2$toker_sex) & (pheno2$text_sex == pheno2$massir_sex), ]
        }
      }
      
      # check that there are enough labels
      filt_pheno <- filt_pheno[filt_pheno$text_sex %in% c("female", "male"),]
      sex_counts <- table(filt_pheno$text_sex)
      print(sex_counts)
      if (sex_counts[['female']] >= NUM_COUNTS & sex_counts[['male']] >= NUM_COUNTS){
        # KEEP --> create class labels
        
        # do we want to discard samples not in this dataset? yes...
        gse.obj2 <- gse.obj 
        gse.obj2$expr <- gse.obj$expr[,filt_pheno$gsm]
        rownames(filt_pheno) <- filt_pheno$gsm
        gse.obj2$class <- as.vector(ifelse(filt_pheno$text_sex=="female", 0,1))
        names(gse.obj2$class) <- filt_pheno$gsm
        gse.obj2$pheno <- filt_pheno
        
        return(gse.obj2)
        
      } else {
        print(sprintf("%s Error - insufficient samples with agreement", gse.id))
        cat(sprintf("%s after sex labeling, insufficient samples with agreement\n", gse.id), file=logfile, append=TRUE)
        return(NA)
      }
    }
  }, error = function(err){ # catch loading errors
    print(sprintf("%s error loading", gse.id))
    cat(sprintf("%s unknown error during sex labeling\n", gse.id), file=logfile, append=TRUE)
    print(err)
    return(NA)
  })
}

runSexLab <- function(fname){
  load.Rdata( sprintf("%s/%s", in_dir, fname), "chunk_ds")
  gse.obj2 <- lapply(chunk_ds, function(x) studySexLabel(x))
  chunk_ds <- gse.obj2
  save(chunk_ds, sprintf("%s/%s", in_dir, fname)) # save the data in the same file
}

args <- commandArgs(trailingOnly=TRUE)
chunk_num <- as.numeric(args[1])
in_dir <- "gses/rObj"
#list.items <- list.files(in_dir)
list.items <- list.files(in_dir, sprintf("chunk%s_", chunk_num))
#out_dir <- "gses/check_sex_lab"
logfile <- sprintf("logs/%s_%s.log", "chunk_label", chunk_num)

for (i in 1:length(list.items)){
  runSexLab(list.items[[i]])
} 

#tmp <- chunkRun(list.items, runSexLab, out_dir, chunk_num, CHUNK_SIZE=100, GROUP_SIZE=10, log_prefix="chunk_label")
  
