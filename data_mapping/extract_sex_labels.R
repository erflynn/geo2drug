# extract_sex_labels.R
# E Flynn
#
# Code for extracting expression-based sex, tissue labels from the drug_mesh data
# run server-side after all of the GSEs have been downloaded as MetaIntegrator objects into a directory

# list the GSE files in the directory
DATA.DIR <- "drug_mesh" #"hc_drug"
list.gse.files <- list.files(DATA.DIR)
list.gses <- sapply(list.gse.files, function(x) strsplit(x, ".", fixed=TRUE)[[1]][1]) # 3184 out of 4706... not great (list: gse_w_drug_mesh.txt)

# function for tissue states 
getMaxTissue <- function(ts, tiss_names){
  apply(ts[,tiss_names], 1, function(y) tiss_names[unlist(which.max(y))])
}

# load the data and grab sex labels
labelMat <- lapply(list.gses, function(gse){
  load(sprintf("%s/%s.RData", DATA.DIR, gse)) # dataObj
  tryCatch({
    ts <- getMaxTissue(dataObj$tissue, colnames(dataObj$tissue)[1:(ncol(dataObj$tissue)-3)])
    df <-rbind(do.call(rbind,c(dataObj$sexlab)), t(ts))  
    rownames(df)[nrow(df)]<-"tissue"
    
  }, error = function(error){
    print(gse)
    ts <- t(rep(NA, length(dataObj$gsms)))
    df <-rbind(do.call(rbind,c(dataObj$sexlab)))
  })
  return(df)
})
names(labelMat) <- list.gses # label with the GSEs

labelMatDF <- do.call(cbind, labelMat)
save(labelMat, file="labelMat.RData")
write.csv(t(labelMatDF), file="expr_label_mat.csv", row.names=TRUE)
