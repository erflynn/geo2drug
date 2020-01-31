require('tidyverse')
require('MetaIntegrator')

MIN.ROWS <- 10000
MIN.CASES <- 10
args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
run_v <- args[2]

all_dat <- read_csv(sprintf("data/01_sample_lists/%s_training.csv", organism))
if (run_v=="full"){
   all_dat2 <- read_csv(sprintf("data/01_sample_lists/%s_training_full.csv", organism))
   all_dat <- rbind(all_dat, all_dat2)
} 
print(nrow(all_dat))

# iterate through the gses
list.gses <- unique(all_dat$gse)
list.train.obj <- lapply(list.gses, 
       function(gse){
         print(gse)
 	 gse_dat <- all_dat[all_dat$gse==gse,]
         gse2 <- strsplit(gse, "-")[[1]][[1]]
         miceadds::load.Rdata(sprintf("data/03_silver_std/%s/00_mat_files/%s_mat.RData", organism, gse2), "mat_obj")
         mat_obj2 <- mat_obj[[sprintf("%s_series_matrix.txt.gz", gse)]]
         gene_df <- apply(mat_obj2$gene_mat[,gse_dat$gsm], c(1,2), as.numeric)
         pheno_df <- mat_obj2$pheno[gse_dat$gsm,]
         class <- ifelse(gse_dat$text_sex=="female", 0, 1)
         names(class) <- gse_dat$gsm
         keys <- rownames(gene_df)
         names(keys) <- keys
         
         # create a meta-object
         my.obj <- list("expr"=gene_df, "pheno"=pheno_df, "class"=class, "keys"=keys, "formattedName"=gse)
	 if (nrow(gene_df) < MIN.ROWS | sum(class==0) < MIN.CASES | sum(class==1) < MIN.CASES ){
	    print(sprintf("%s had too few rows, n=%s.", gse, nrow(gene_df)))
	    print(table(class))
	    return(NA)
	 }
         if(!checkDataObject(my.obj, "Dataset")){
           print(sprintf("Error with %s", gse))
           return(NA)
         } 

         #obj <- list(gse=my.obj)
         return(my.obj)
})

names(list.train.obj) <- list.gses

# now run meta-integrator on this!
metaObject <- list("originalData"=list.train.obj[!is.na(list.train.obj)])

save(metaObject, file=sprintf("data/03_silver_std/%s/04_meta_res/input_metaObj_%s.RData", organism, run_v))
metaObject <- runMetaAnalysis(metaObject)
save(metaObject, file=sprintf("data/03_silver_std/%s/04_meta_res/metaObj_%s.RData", organism, run_v)  )


