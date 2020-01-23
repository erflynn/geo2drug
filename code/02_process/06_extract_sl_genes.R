
# the goal of this script is to extract the sex labeling genes from the meta-object
# // TODO - selection of filter effects results, examine this


require('MetaIntegrator')
require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
run_v <- args[2]

miceadds::load.Rdata(sprintf("data/03_silver_std/%s/04_meta_res/metaObj_%s.RData", organism, run_v), "metaObj") 

num_studies <- length(metaObj$originalData)

if (run_v == "full"){ # remove loo b/c so big
  metaCopy <- metaObj
  loo <- metaCopy$leaveOneOutAnalysis
  loo2 <- loo[sapply(loo, function(x) length(x))!=0]
  metaCopy$leaveOneOutAnalysis <- loo2
  metaObj <- metaCopy
}

## try a couple different sets of filters for each so we have all these data
# diff filters for each to get a reasonable number of genes
if (organism == "rat"){
  metaObj <- filterGenes(metaObj, isLeaveOneOut=TRUE, 
                         FDRThresh = 0.05, effectSizeThresh=0.5)
}
if (organism=="mouse"){
  metaObj <- filterGenes(metaObj, isLeaveOneOut=FALSE, 
                         FDRThresh = 0.05, effectSizeThresh=0.4,
                         numberStudiesThresh = 10)
  
} else {
  metaObj <- filterGenes(metaObj, isLeaveOneOut=FALSE, 
                         FDRThresh = 0.05, effectSizeThresh=0.5,
                         numberStudiesThresh = floor(0.6*num_studies))
}

res <- metaObj$filterResults[[1]]

m_genes <- res$posGeneNames # MALES
f_genes <- res$negGeneNames # FEMALES


# grab chromosome name + other symbols
miceadds::load.Rdata(sprintf("gpl_ref/%s_gene_map.RData", organism), "gene_map")
if (organism != "human"){
  sel.list <- c( "entrezgene_id", "chromosome_name")
} else {
  sel.list <- c("hgnc_symbol", "entrezgene_id", "chromosome_name")
}
mgenes_df <- gene_map %>% 
  filter(entrezgene_id %in% m_genes) %>% 
  select(sel.list) %>% 
  unique() %>%
  filter(chromosome_name %in% c(1:22, "X", "Y"))

fgenes_df <- gene_map %>% 
  filter(entrezgene_id %in% f_genes) %>% 
  select(sel.list) %>% 
  unique() %>% 
  filter(chromosome_name %in% c(1:22, "X", "Y"))

fgenes_df %>% write_csv(sprintf("data/03_silver_std/%s/04_meta_res/fgenes_%s.csv", organism, run_v))
mgenes_df %>% write_csv(sprintf("data/03_silver_std/%s/04_meta_res/mgenes_%s.csv", organism, run_v))

# write out the FDR info
res_summary <- summarizeFilterResults(metaObj, getMostRecentFilter(metaObj))
res_summary$pos %>% write_csv(sprintf("data/03_silver_std/%s/04_meta_res/mgenes_summary_%s.csv", organism, run_v))
res_summary$neg %>% write_csv(sprintf("data/03_silver_std/%s/04_meta_res/fgenes_summary_%s.csv", organism, run_v))
