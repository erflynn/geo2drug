
# load some of the recount data


## Load library
library('recount')

## Find a project of interest
project_info <- abstract_search('GSE32465')

## Download the gene-level RangedSummarizedExperiment data
download_study(project_info$project)

load(file.path(project_info$project, 'rse_gene.Rdata'))
colData(rse_gene)$geo_accession
dim(assays(rse_gene)@listData$counts)

head(assays(rse_gene)@listData$counts) # <-- this is the expression matrix

rse_gene$characteristics


rse_gene$title
rse_gene$characteristics

# http://metasra.biostat.wisc.edu/api/v01/samples.json?study=SRP009615
# --> this gives the metadata
# wou