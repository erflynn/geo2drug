# check_gpl_success.R
# E Flynn
# 10/1/2019
# 
# Code for examining whether the GPL succeeds at downloading and the output looks good


# load in the GPL counts list
hs_count_gpls <- read_csv("data/sample_lists/count_gpl_human.csv")
hs_count2 <- filter(hs_count_gpls, n > 1)

head(hs_count2$gpl)



# list the GPLs in the directory
#  ls -l > ../list_human_downloaded.txt
hs_download_dump <- read_table("data/sample_lists/list_human_downloaded.txt", skip=1, col_names=FALSE)

gpl_downloaded <- sapply(hs_download_dump$X9, function(x) strsplit(x, ".", fixed=TRUE)[[1]][[1]])

hs_plus_download <- hs_count2 %>%
  mutate(downloaded=(gpl %in% gpl_downloaded))

# this is not promising...
ggplot(hs_plus_download, aes(y=log(n), x=as.factor(downloaded)))+geom_violin()+geom_point()

sum(hs_plus_download[hs_plus_download$downloaded,]$n) # 10329
sum(hs_plus_download[!hs_plus_download$downloaded,]$n) # 8178
filter(hs_plus_download, n > 100) %>% group_by(downloaded) %>% count()

# plot something - do we care?