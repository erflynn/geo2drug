require('tidyverse')
require('MetaIntegrator')



all_dat <- read_csv("../data/sample_lists/human_training.csv")

# iterate through the gses
list.gses <- unique(all_dat$gse)
list.train.obj <- lapply(list.gses, 
       function(gse){
         print(gse)
 gse_dat <- all_dat[all_dat$gse==gse,]
         gse2 <- strsplit(gse, "-")[[1]][[1]]
         miceadds::load.Rdata(sprintf("../data/silver_std/%s_mat.RData", gse2), "mat_obj")
         mat_obj2 <- mat_obj[[sprintf("%s_series_matrix.txt.gz", gse)]]
         gene_df <- apply(mat_obj2$gene_mat[,gse_dat$gsm], c(1,2), as.numeric)
         pheno_df <- mat_obj2$pheno[gse_dat$gsm,]
         class <- ifelse(gse_dat$text_sex=="female", 0, 1)
         names(class) <- gse_dat$gsm
         keys <- rownames(gene_df)
         names(keys) <- keys
         
         # create a meta-object
         my.obj <- list("expr"=gene_df, "pheno"=pheno_df, "class"=class, "keys"=keys, "formattedName"=gse)
         if(!checkDataObject(my.obj, "Dataset")){
           print(sprintf("Error with %s", gse))
           return(NA)
         }
         #obj <- list(gse=my.obj)
         return(my.obj)
})

names(list.train.obj) <- list.gses

# now run meta-integrator on this!
metaObject <- list("originalData"=list.train.obj[!is.na(list.train.obj))

#load("../data/input_metaObj.RData")
metaObject <- runMetaAnalysis(metaObject)
save(metaObject, file="../data/metaObj.RData")  


