
# to run:
# Rscript gse_to_mat.R gse.file out.dir idx

require('exprsex')
require('MetaIntegrator')

GSE.DIR <- "/scratch/users/erflynn/sex_labeling/geo_pipeline/gses/human/matrix/"
GPL.DIR <- "/scratch/users/erflynn/sex_labeling/geo_pipeline/gpl_ref/human/"
SIZE.CHUNK <- 50

# parse arguments
args <- commandArgs(trailingOnly=TRUE)
gse.file <- args[1]
gse.list <- read.csv(gse.file, header=TRUE)
OUT.DIR <- args[2]
idx <- as.numeric(args[3])
print(idx)

NUM.CHUNKS <- ceiling(nrow(gse.list)/SIZE.CHUNK)
end_idx <- ifelse((NUM.CHUNKS-1) == idx ,nrow(gse.list), (idx+1)*SIZE.CHUNK)
gse.list <- gse.list[(idx*SIZE.CHUNK):end_idx,]

print(gse.list[,1])

gses.to.run <- setdiff(gse.list[,1], c("GSE39144", "GSE76246", "GSE79945", "GSE76516", "GSE70564", "GSE25219","GSE37138", "GSE18927", "GSE30727", "GSE50421", "GSE40492", "GSE31983","GSE19090", "GSE26106","GSE70565", "GSE76519", "GSE84890" , "GSE77714", "GSE28387"))

lapply(gses.to.run, function(gse){

tryCatch({
	print(gse)
	 

	 out.file <- sprintf("%s/%s_mat.RData", OUT.DIR, gse)
	 if (!file.exists(out.file)){
	    geo.obj <- exprsex::getPrepGSE(gse, gpl.dir=GPL.DIR, gse.dir=GSE.DIR)
	    save(geo.obj, file=out.file)
   	 } else {
	 print("already loaded")
	 }

}, err= function(e){
   print(sprintf("%s error loading", gse))
   print(e)
})
})