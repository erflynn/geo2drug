# E Flynn
# Sex breakdown of *metadata annotations* (using ALE) for ALL mouse and rat studies in GEO

require('tidyverse')

reformatAleDat <- function(dat){
  dat2 <- dat %>% rename("gsm"="X", "gse"="ExperimentID", "gpl"="PlatformID") %>% 
    mutate(gsm=paste("GSM", gsm, sep=""), gse=paste("GSE", gse, sep=""), gpl=paste("GPL", gpl, sep=""))
  return(dat2)
}
summarizeToStudy <- function(dat){
  study_dat <- dat %>% group_by(gse) %>% 
    summarise(num_f=sum(Gender=="F"), num_m=sum(Gender=="M"), total=n(), gpl=paste(unique(gpl), collapse="|")) %>% 
    mutate(study_type = ifelse(num_f == 0, ifelse(num_m ==0, "no_labels", "m"), ifelse(num_m==0, "f", "both")))
  return(study_dat)
}

args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
ale.file <- args[2]

dat <- read.csv(ale.file, stringsAsFactors=FALSE)
#mouse_dat <- read.csv(sprintf("%s/ale_out_10090", ALE.DIR), stringsAsFactors = FALSE)
#rat_dat <- read.csv(sprintf("%s/ale_out_10116", ALE.DIR), stringsAsFactors = FALSE)
#human_dat <- read.csv(sprintf("%s/ale_out", ALE.DIR), stringsAsFactors = FALSE)

# reformat the data to get the by-sample breakdown
dat2 <- reformatAleDat(dat)
print(table(dat2$Gender)  )

# write out the sex labels for use later
dat3 <- dat2 %>% select(gse, gsm, gpl, Gender) %>% filter(Gender!="")
length(unique(dat3$gse))
ale_combined <- read_csv("data/deprecated/ale_combined_data.csv")
ale2 <- ale_combined %>% filter(text_sex !="")
ale_combined[ale_combined$gse=="GSE112371" & !is.na(ale_combined$gse),]

dat3 %>% write_csv(sprintf("data/01_sample_lists/%s_ale_sex_lab.csv", organism))

# summarize to get the by-study breakdown
dat_study <- summarizeToStudy(dat2)

print(table(dat_study$study_type))

write.csv(dat_study, file=sprintf("data/01_sample_lists/%s_gse_counts.csv", organism), row.names=FALSE)
