# methods for de-duplicating a list of GSE/GSM


isSuper <- function(s1.samples, s2.samples){
  # returns: 0 (same), 1, 2, or NA
  if (union(s1.samples, s2.samples)==intersect(s1.samples, s2.samples)){
    return(0)
  }
  
  if (length(setdiff(s2.samples, s1.samples)) == 0){
    return(2)
  }
  if (length(setdiff(s1.samples, s2.samples)) == 0){
    return(1)
  } 
  return(NA)
}


getEarlierStudy <- function(id1, id2, study_info){
  dates <- sapply(c(id1, id2), function(x) 
    study_info$submission_date[study_info$gse==x])
  earliest.study <- c(id1, id2)[which.min(dates)]
}

deduplicate <- function(input_file, out_file){
  my_studies <- read.csv(input_file, stringsAsFactors = FALSE)
  my_studies <- my_studies[!duplicated(my_studies[,c("gse", "gsm")]),]
  gsm.to.gse <- split(my_studies$gse, my_studies$gsm)
  table(sapply(gsm.to.gse, length)) # many GSMs are in multiple!
  mult.studies <- gsm.to.gse[sapply(gsm.to.gse, length) > 1]
  mult.study.strs <- lapply(mult.studies, function(x) paste(unique(x[order(x)]), collapse=" "))
  study.comb <- unique(mult.study.strs) 

  gse.to.gsm <- split(my_studies$gsm, my_studies$gse)
  list_gses_form <- paste(unique(my_studies$gse), collapse="\', \'")
  require('GEOmetadb')
  con <- dbConnect(SQLite(),'../GEOmetadb.sqlite') 
  study_info <- dbGetQuery(con, sprintf("SELECT gse, submission_date FROM gse WHERE gse IN ('%s');", list_gses_form))
  dbDisconnect(con)
  study_info$submission_date <- as.Date(study_info$submission_date)
  
  removed_gses <- c()
  removed_gsms <- list()
  for (comb in study.comb){
    list_comb <- strsplit(comb, " ", fixed=TRUE)[[1]]
    for (i in 1:(length(list_comb)-1)){
      if (list_comb[i] %in% removed_gses){
        print(sprintf("%s already removed", list_comb[i]))
        break 
      }
      for (j in (i+1):length(list_comb)){
        pair <- c(list_comb[i], list_comb[j])
        print(pair)
        id1 <- pair[[1]]
        id2 <- pair[[2]]
        if (id2 %in% removed_gses){
          print(sprintf("%s already removed", id2))
          break 
        }
        s1.samples <- gse.to.gsm[[id1]]
        s2.samples <- gse.to.gsm[[id2]]
        super <- isSuper(s1.samples, s2.samples) 
        if (is.na(super)){
          id.early <- getEarlierStudy(id1, id2, study_info)
          id.other <- pair[pair!=id.early]
          overlap <- intersect(s1.samples, s2.samples)
          print(sprintf("The overlap between %s and %s is %s, assigning to %s", id1, id2, length(overlap), id.early))
          
          removed_gsms[[sprintf("%s", id.other)]] <- list(overlap)
        } else {
          # keep one
          if (super == 0){
            id.early <- getEarlierStudy(id1, id2, study_info)
            id.other <- pair[pair!=id.early]
            print(sprintf(" %s and %s are the same, keeping %s", id1, id2, id.early))
            removed_gses <- c(removed_gses, id.other)
          }
          # keep the superset
          else {
            idsup <- pair[[super]]
            id.other <- pair[[-super]]
            print(sprintf(" %s is a superset of %s", idsup, id.other))
            removed_gses <- c(removed_gses, id.other)
          }
        }
      }
    }
  }
  
  my.studies2 <- my_studies[!my_studies$gse %in% removed_gses,]
  my.studies3 <- my.studies2
  for (i in 1:length(removed_gsms)){
    gse <- names(removed_gsms)[i]
    gsms <- removed_gsms[[i]][[1]]
    my.studies3 <- my.studies3[!(my.studies3$gse ==gse & my.studies3$gsm %in% gsms),]
  }
  
  write.csv(my.studies3, file=out_file, row.names=FALSE)
}