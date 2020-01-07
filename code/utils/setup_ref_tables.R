
require('tidyverse')

setup_ref_tables <- function(){
  require('biomaRt')
  ensembl_h <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl") 
  biomaRt::listAttributes(ensembl_h) %>% 
    dplyr::select(name) %>% 
    dplyr::filter(str_detect(name, "uni"))
  human.dat.map <- biomaRt::getBM(attributes= c("ensembl_gene_id", "hgnc_symbol", "refseq_mrna", "entrezgene_id", "chromosome_name"), mart=ensembl_h) 
  dat.map <- human.dat.map
  save(dat.map, file="../../data/human_gene_map.RData")
  
  ensembl_m <- useMart("ensembl", dataset="mmusculus_gene_ensembl") 
  mouse.dat.map <- getBM(attributes= c("ensembl_gene_id", "hgnc_symbol", "refseq_mrna", "entrezgene_id", "chromosome_name"), mart=ensembl_m) 
  dat.map <- mouse.dat.map
  save(dat.map, file="../../data/mouse_gene_map.RData")
  
  ensembl_r <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl") 
  rat.dat.map <- getBM(attributes= c("ensembl_gene_id", "refseq_mrna", "entrezgene_id", "chromosome_name"), mart=ensembl_r) 
  dat.map <- rat.dat.map
  save(dat.map, file="../../data/rat_gene_map.RData")
}

extract_entrez <- function(organism){
  load(sprintf("../../data/%s_gene_map.RData", organism))
  entrez_df <- dat.map %>% dplyr::filter(entrezgene_id != "") %>% 
    dplyr::select(entrezgene_id) %>% unique()
  entrezids <- sapply(entrez_df$entrezgene_id, as.character)
  save(entrezids, file=sprintf("../../data/%s_entrezids.RData", organism))
}


get_unigene <- function(organism){
  require('mygene')
  
  load(sprintf("data/%s_entrezids.RData", organism))
  df_unigene <- queryMany(entrezids, scopes = "entrezgene", fields = "unigene", species = organism, returnall = TRUE)
  
  table(sapply(df_unigene$response[,"unigene"], is.null)) 
  df_unigene$response$unigene <- sapply(df_unigene$response$unigene, function(x) paste(x, collapse=";"))
  my_df <- df_unigene$response[df_unigene$response$unigene!="",]
  unigene <- data.frame(my_df) %>% 
    dplyr::select("query", "unigene") %>% 
    dplyr::rename(entrezgene_id=query) %>%
    tidyr::separate_rows(unigene, sep=";")
  print(head(unigene))
  save(unigene,
       file=sprintf("data/%s_unigene.RData", organism ))
}


get_genbank <- function(organism){
  
  if (organism=="human"){
    library(org.Hs.eg.db)
    x <- org.Hs.egACCNUM
  } else if (organism == "mouse"){
    require("org.Mm.eg.db")
    x <- org.Mm.egACCNUM
  } else if (organism == "rat"){
    require("org.Rn.eg.db")
    x <- org.Rn.egACCNUM
  } else {
    print("Only implemented for human, rat, and mouse")
    return(NA)
  }
  
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  collapsed <- lapply(xx, function(y) paste(unique(y), collapse=" /// "))
  df <- data.frame(cbind("entrezgene_id"=names(xx), "genbank"=collapsed)  )
  rownames(df) <- NULL
  genbank <- df %>% separate_rows(genbank, sep=" /// ")
  save(genbank, file=sprintf("data/%s_genbank.RData", organism))
}





