

getSexLabStats <- function(df){
  # function to get stats of numbers m,f for a particular study 
  count.per.study <- df %>% group_by(gse) %>% summarize(num_f=sum(sexlab=="female"), num_m=sum(sexlab=="male"))
  m.only <- filter(count.per.study, num_f==0)
  f.only <- filter(count.per.study, num_m==0)
  mixed.sex <- filter(count.per.study, num_m!=0, num_f!=0) 
  
  num_studies <- length(unique(comb_labels$gse))
  list.stats <- list("f"=nrow(f.only), "m"=nrow(m.only), "both"=nrow(mixed.sex), "unlab"=num_studies-length(unique(df$gse)))
  return(list.stats)
}


vec_to_df <- function(vec, item_sep="|") {
  # converts a named vector to a data frame
  # TODO - clean this up and make it a better function, speed up!
  #   - column names are HIDEOUS
  
  df <- data.frame(do.call(rbind, lapply(vec, function(x) paste(x, collapse=item_sep)) ))
  df$id <- names(vec)
  return(df)
}

grabMeshType <-  function(x) {substr(x, 1, 1)} 
# gets the first character