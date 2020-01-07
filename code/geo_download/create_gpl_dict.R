# create_gpl_dict.R
#

require('MetaIntegrator')
require('tidyverse')
gse.name <- "GSE10846"
gse.obj <- getGEOData(gse.name)$originalData[[1]]

# probes to genes conversion
keys <- gse.obj$keys
#list.keys <- keys[keys %in% gene.list]
list.keys <- keys
key.df <- data.frame(list.keys, names(list.keys))
colnames(key.df) <- c("gene", "probes")
key.df$gene <- as.character(key.df$gene)
key.df$probes <- as.character(key.df$probes)
key.df2 <- separate_rows(key.df, gene, sep=",")

gene.to.probe <- split(key.df2$probes,  key.df2$gene)

