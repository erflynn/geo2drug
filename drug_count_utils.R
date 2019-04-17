

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


# vec_to_df <- function(vec, item_sep="|") {
#   # converts a named vector to a data frame
#   # TODO - clean this up and make it a better function, speed up!
#   #   - column names are HIDEOUS
#   
#   df <- data.frame(do.call(rbind, lapply(vec, function(x) paste(x, collapse=item_sep)) ))
#   df$id <- names(vec)
#   return(df)
# }

vec_to_df <- function(vec_info, item_sep="|"){
  df <- do.call(rbind, 
          lapply(vec_info, function(x) 
            data.frame(lapply(x, function(y) paste(y, collapse=item_sep)), stringsAsFactors=FALSE)))
  df$id <- names(vec_info)
  return(df)
}

grabMeshType <-  function(x) {substr(x, 1, 1)} 
# gets the first character

createNameSynVocab <- function(df, id_col){
  # Creates a controlled vocabulary from the name and synonyms, mapping them to the id_col
  # Assumes the df has the columns: name, synonyms, <id_col>
  # output is a df with format name, <id_col>
  name_syn <- rbind(df[,c("name", id_col)], rename(df[,c("synonyms", id_col)], "name"="synonyms"))
  name_syn <- separate_rows(name_syn, name, sep="\\|")
  
  # convert to lower
  name_syn$name <- sapply(name_syn$name, function(x) tolower(str_trim(x)))
  name_syn <- name_syn[!duplicated(name_syn),]
  
  # remove empty
  name_syn <- filter(name_syn, id_col!="")
  name_syn <- filter(name_syn, name!="")
  
  return(name_syn)
}



labelDrugGSE <- function(gse_metadata){
  # returns a drugbank ID of the drug in the gse metadata, or NA if none is present
  
}

labelCellGSE <- function(gse_metadata){
  # returns the cell line name of the data
}

