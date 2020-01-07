

ORG.ENSEMBL.SUFF <- list("human"="G", "mouse"="MUSG", "rat"="RNOG")
ORG.UNIGENE.PRE <- list("human"="Hs", "mouse"="Mm", "rat"="Rn")

reform_entrez_df <- function(df, my.col.name="gene"){
  df[,my.col.name] <- sapply(df[,my.col.name], function(x) 
    paste(stringr::str_extract_all(x, "[0-9]+")[[1]], collapse=";"))
  df %>% tidyr::separate_rows(my.col.name, sep=";")
}


find_col_loc <- function(df, my.col, ex.str){
  ex_rows <- which(stringr::str_detect(df[,my.col], ex.str))
  my_str <- df[ex_rows[[1]], my.col]
   
  mult_fields <- stringr::str_split(my_str, " /// ")[[1]]
  lst_fields <- stringr::str_split(mult_fields, " // ")
  idces <- lapply(lst_fields, function(x)
    which(sapply(x, function(y)
      stringr::str_detect(y, sprintf("^%s", ex.str)))))
  idx <- as.numeric(unique(unlist(idces))[[1]])
  
  gene_vals <- lapply(df[,my.col], function(x){
    mult_fields <- stringr::str_split(x, " /// ")[[1]]
    lst_fields <- stringr::str_split(mult_fields, " // ")
    genes <- unique(lapply(lst_fields, function(x) x[[idx]]))
    genes <- genes[genes!="---"]
    paste(genes, collapse=" /// ")
  })
  df2 <- data.frame(cbind("probe"=df[,1], "gene_col"=gene_vals))
  rownames(df2) <- NULL
  df3 <- df2 %>% tidyr::separate_rows(gene_col, sep=" /// ")
  return(df3)
}

parse_multi_col <- function(df, my.col, pattern){
  # check whether it matches the whole pattern
  if (any(stringr::str_detect(df[,my.col], sprintf("^%s$", pattern)))){
    df2 <- df[,c(1, my.col)]
    colnames(df2) <- c("probe", "gene_col")
    df3 <- df2 %>% tidyr::separate_rows(gene, sep=" /// ")
    return(df3)
  }
  gene_vals <- str_extract_all(df[,my.col], pattern)
  gene_vals_col <- lapply(gene_vals, function(x) 
    paste(x, collapse=" /// "))
  probe.ids <- df[,1]
  df2 <- data.frame(cbind("probe"=probe.ids, "gene_col"=gene_vals_col))
  df3 <- df2 %>% tidyr::separate_rows(gene_col, sep=" /// ")
  return(df3)
}

check_entrez_overlap <- function(gpl.df, cols, organism){
  # grab the entrezids associated with the organism
  load(sprintf("data/%s_entrezids.RData", organism))
  
  overlap.lengths <- sapply(1:length(cols), function(i){
    length(intersect(entrezids, stringr::str_extract_all(gpl.df[,cols[i]], "[0-9]+")))
  })
  names(overlap.lengths) <- cols
  return(overlap.lengths)
}

map_from_refseq <- function(gpl.df, my.col, organism){
  # check that it matches the whole column -- if not, parse
  probe.gene <- find_col_loc(gpl.df, my.col, "N[R|M][_][\\d]+[_.-]*[\\w\\d]*")
  colnames(probe.gene) <- c("probe", "refseq_mrna")
  load(sprintf("data/%s_gene_map.RData", organism))  
  ref_entr <- dat.map %>% 
    dplyr::select(refseq_mrna, entrezgene_id) %>% 
    dplyr::filter(!is.na(entrezgene_id)) %>%
    unique()
  dplyr::inner_join(probe.gene, ref_entr) %>% 
    dplyr::rename(gene=entrezgene_id)
}

map_from_ensembl <- function(gpl.df, my.col, organism){
  probe.gene <- parse_multi_col(gpl.df, my.col, sprintf("[%s][\\d]+[.-]*[\\w\\d]*", org.ensembl.id))
  colnames(probe.gene) <- c("probe", "ensembl_gene_id")
  load(sprintf("data/%s_gene_map.RData", organism))
  ens_entr <- dat.map %>% 
    dplyr::select(ensembl_gene_id, entrezgene_id) %>% 
    dplyr::filter(!is.na(entrezgene_id)) %>%
    unique()
  dplyr::inner_join(probe.gene, ens_entr) %>% dplyr::rename(gene=entrezgene_id)
}

map_from_unigene <- function(gpl.df, my.col, organism){
  probe.gene <- parse_multi_col(gpl.df, my.col, sprintf("%s.[0-9]+",ORG.UNIGENE.PRE[[organism]]))
  colnames(probe.gene) <- c("probe", "unigene")
  load(sprintf("data/%s_unigene.RData", organism))
  dplyr::inner_join(probe.gene, unigene) %>% dplyr::rename(gene=entrezgene_id)
}


map_from_hgnc <- function(gpl.df, my.col, organism){
  if (organism=="rat"){
    print("error no HGNC symbol mapping for rat")
    return(NA)
  }
  probe.gene <- find_col_loc(gpl.df, my.col, "GAPDH|gapdh|Gapdh")
  colnames(probe.gene) <- c("probe", "hgnc_symbol")
  load(sprintf("data/%s_gene_map.RData", organism))  
  hgnc_entr <- dat.map %>% 
    dplyr::select(hgnc_symbol, entrezgene_id) %>% 
    dplyr::filter(!is.na(entrezgene_id)) %>%
    unique()
  dplyr::inner_join(probe.gene, ref_entr) %>% dplyr::rename(gene=entrezgene_id)
}

map_from_genbank <- function(gpl.df, my.col, organism, genbank_str){
  probe.gene <- find_col_loc(gpl.df, my.col, genbank_str)
  colnames(probe.gene) <- c("probe", "genbank")
  load(sprintf("data/%s_genbank.RData", organism))
  dplyr::inner_join(probe.gene, genbank) %>% dplyr::rename(gene=entrezgene_id)
}


parse_entrez_from_gpl <- function(gpl.name, MIN.OVERLAP=8000){

  # read in GPL data
  gpl <- GEOquery::getGEO(gpl.name)

  # check the organism
  organism <- gpl@header$taxid
  if (stringr::str_detect(organism, "9606")){
    org.name <- "human"
  } else if (stringr::str_detect(organism, "10090")){ 
    org.name <- "mouse"
  } else if (stringr::str_detect(organism, "10116")){
    org.name <- "rat"    
  } else {
    print(sprintf("Error for %s, gpl key extraction is only implemented for human, rat, and mouse.", 
                  gpl_name))
    return(NA)
  }

  # TODO - what if the probe column is the row name
  gpl.df <- gpl@dataTable@table
  probe.ids <- gpl.df[,1]
  
  # ------ Find columns labeled ENTREZ IDs ------ #
  entrez.col1 <- which(grepl("entrez", 
                             gpl@dataTable@columns$Column, ignore.case = TRUE))
  entrez.col2 <- which(grepl("entrez", 
                             gpl@dataTable@columns$Description, ignore.case = TRUE))
  entrez.col <- unique(entrez.col1, entrez.col2)
  if (length(entrez.col) != 0){
    # look for overlap
    overlap.lengths <- check_entrez_overlap(gpl.df, entrez.col, org.name)
    if (max(overlap.lengths) > MIN.OVERLAP){
      entrez.col <- as.numeric(names(overlap.lengths)[which.max(overlap.lengths)[[1]]])
      df <- data.frame("probe"=probe.ids, "gene"=gpl.df[,entrez.col])
      return(reform_entrez_df(df))
    }
  }
  
  # ------ Find any all integer columns, check for overlap ------ #
  int_only_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "^[0-9]+$"), na.rm=TRUE)
  })
  int_only_cols <- (int_only_cols > MIN.OVERLAP)
  col_idx <- which(int_only_cols)
  if (length(col_idx) != 0){
    overlap.lengths <- check_entrez_overlap(gpl.df, col_idx, org.name)
    if (max(overlap.lengths) > MIN.OVERLAP){
      entrez.col <- as.numeric(names(overlap.lengths)[which.max(overlap.lengths)[[1]]])
      df <- data.frame("probe"=probe.ids, "gene"=gpl.df[,entrez.col])
      return(reform_entrez_df(df))
    }
  }
  
  # ------ MAP FROM REFSEQ ------ #
  refseq_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "NM_[0-9]+"), na.rm=TRUE)
  })
  if (max(refseq_cols) > MIN.OVERLAP){
    refseq_col <- which.max(refseq_cols)[[1]]
    print("parsed from refseq")
    return(map_from_refseq(gpl.df, refseq_col, org.name))
  }
  
  # ------ MAP FROM ENSEMBL ------ #
  ensembl_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], sprintf("ENS%s[0-9]+", ORG.ENSEMBL.SUFF[["organism"]])), na.rm=TRUE)
  })
  if (max(ensembl_cols) > MIN.OVERLAP){
    ensembl_col <- which.max(ensembl_cols)
    print("parsed from ensembl")
    return(map_from_ensembl(gpl.df, ensembl_col, org.name))
  }
  
  # ------ MAP FROM UNIGENE ----- #
  unigene_cols <- sapply(1:ncol(gpl.df), function(i){
    sum(stringr::str_detect(gpl.df[,i], "Hs.|Mm.|Rn."), na.rm=TRUE)
  })
  if (max(unigene_cols) > MIN.OVERLAP){
    unigene_col <- which.max(unigene_cols)
    print("parsed from unigene")
    return(map_from_unigene(gpl.df, unigene_col, org.name))
  }

  # ------  MAP FROM HGNC  ------ #
  hgnc_cols <- sapply(1:ncol(gpl.df), function(i){
      grepl("GAPDH", gpl.df[,i], ignore.case=TRUE)
  })
  if (length(hgnc_cols) > 0){
    print("parsed from HGNC")
    return(map_from_hgnc(gpl.df, hgnc_cols[[1]], org.name))
    
  }

  # ------ MAP FROM GenBank ------ #
  #rat_genbank <- paste((genbank %>% filter(entrezgene_id=="24383"))[1:5,"genbank"], collapse="|")
  #mouse_genbank <- paste((genbank %>% filter(entrezgene_id=="14433"))[1:5,"genbank"], collapse="|")
  #human_genbank <- paste((genbank %>% filter(entrezgene_id=="2597"))[1:5,"genbank"], collapse="|")
  LIST.GENBANK.STR <- 
    list("rat"="AAA40814|AAA41193|AAA41795|AAB19105|AAH59110",
         "mouse" ="AAA37659|AAH82592|AAH83065|AAH83079|AAH83080",
         "human"="AAA52496|AAA52518|AAA52519|AAA53191|AAA86283")
  genbank_str <- LIST.GENBANK.STR[[org.name]]
  genbank_cols <- sapply(1:ncol(gpl.df), function(i){
    grepl(genbank_str, gpl.df[,i])
   })
  if (length(genbank_cols) > 0){
   print("parsed from GenBank")
    return(map_from_genbank(gpl.df, genbank_cols[[1]], org.name, genbank_str))
  }
  print("No mapping")
  return(NA)
}
parsed2 <- parse_entrez_from_gpl("GPL10558")
parsed3 <- parse_entrez_from_gpl("GPL17586")
