
#' Wrapper for running massiR.
#' 
#' @param expData expression data matrix with rows as genes and columns as samples
#' @param ychr.genes a list of y chromosome genes
#' @param threshold the threshold parameter for massiR, defaults to 3
#' @param na.cut max fraction missingness for the sample data, defaults to 0.3
#' @param min.samples minimum number of samples needed to proceed, defaults to 10
#' 
#' @return a list of sex labels, named by samples
massiRAcc <- function(expData, ychr.genes, threshold=3, na.cut=0.3, min.samples=10){
  # note - this reorders the results!, have fixed this
  
  require('massiR')
  if (ncol(expData)==1 | is.null(colnames(expData))){
    print("Error, massir method - there are not enough samples without missing data")
    return(NULL)
  }
  cols <- colnames(expData)
  expData <- expData[,order(cols)]
  
  # label sex
  ychr.genes2 <- sapply(intersect(ychr.genes, rownames(expData)), as.character)
  
  # remove columns with missing data
  na.count <- apply(expData[ychr.genes2,], 2, function(x) sum(is.na(x)))
  genes.df <- expData[,na.count < na.cut*length(ychr.genes2)]
  
  # fails if there are few columns left
  if (ncol(genes.df) < min.samples){
    print("Error, massir method - there are not enough samples without missing data")
    return(NULL)
  }
  
  cols <- colnames(expData)
  
  # remove genes with missingness
  genes.df2 <- genes.df[complete.cases(genes.df),]
  ychr.genes3 <- sapply(intersect(ychr.genes, rownames(expData)), as.character)
  
  # return an error if few genes are present
  if (length(ychr.genes3)<2){
    print("Error, less than 2 ychr genes present")
    return(NULL)
  }
  ychr.df <- data.frame(row.names=(ychr.genes3))
  ds.y.out <- massi_y(data.frame(genes.df2), ychr.df)
  
  massi.select.out <-
    + massi_select(data.frame(genes.df2), ychr.df, threshold=threshold) 
  
  results <- massi_cluster(massi.select.out)
  sample.results <- data.frame(results$massi.results)
  sex.lab <- sample.results$sex
  names(sex.lab) <- sample.results$ID
  
  # add in the missing columns
  rem.cols <- setdiff( cols, names(sex.lab))
  rem.vec <- rep(NA, length(rem.cols))
  names(rem.vec) <- rem.cols
  sex_lab2 <- c(sex.lab,rem.vec)
  return(sex_lab2[cols])
}


#' Run toker sex labeling using the Toker et al method
#' @param expData expression matrix with rows as genes and columns as samples
#' @param f.genes ids for f-spec genes (XIST)
#' @param m.genes ids for m-spec genes (KDM5D, RPS4Y1)
#' @param min.samples the sample cutoff for returning results, default to 8
#' 
#' @return list of sex labels, NULL if sex labeling failed
tokerSexLab <- function(expData, f.genes="7503", m.genes=c("6192", "8284"),
                        min.samples=8){
  
  keys <- rownames(expData)
  f.genes2 <- sapply(intersect(keys, f.genes), as.character)
  m.genes2 <- sapply(intersect(keys, m.genes), as.character)
  
  if (length(f.genes2)==0 | length(m.genes2)==0){
    print("Error, Toker sex labeling genes not present in dataset")
    return(NULL)
  }
  
  
  genes.df <- expData[c(f.genes2, m.genes2),]
  cols <- colnames(genes.df)
  
  # remove columns with any NAs
  na.cols <- apply(genes.df, 2, function(x) any(is.na(x)))
  genes.df2 <- genes.df[,!na.cols]
  
  # fails if there are few columns left
  if (ncol(genes.df2) < min.samples){
    print("Error, Toker method - there are not enough samples without missing data")
    return(NULL)
  }
  
  clus <- kmeans(t(genes.df2), centers=2)
  sex_lab <- clus$cluster

  # get clusters
  if (length(m.genes2)==1){
    m.clus <- unlist(which.max(clus$centers[,c(m.genes2)])) 
  } else {
    m.clus <- unlist(which.max(rowMeans(clus$centers[,c(m.genes2)])))
  }
  
  if (length(f.genes2)==1){
    f.clus <- unlist(which.max(clus$centers[,c(f.genes2)])) 
  } else {
    f.clus <- unlist(which.max(rowMeans(clus$centers[,c(f.genes2)])))
  }
  
  if (f.clus==m.clus){
    print("Error, Toker method - the same cluster has higher f and m values so we cannot assign.")
    return(NULL)
  }
  sex_lab[sex_lab==f.clus] <- "female"
  sex_lab[sex_lab==m.clus] <- "male"
  
  # add in the missing columns
  rem.cols <- setdiff( cols, names(sex_lab))
  rem.vec <- rep(NA, length(rem.cols))
  names(rem.vec) <- rem.cols
  sex_lab2 <- c(sex_lab,rem.vec)
  return(sex_lab2[cols])
}
