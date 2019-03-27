# condenseMeSH.R
# E Flynn
# 2/27/2019
# 
# Code for condensing MeSH data into a single dataframe

require('rjson')
require('tidyverse')
source('drug_count_utils.R')

mesh_info <- fromJSON(file="data/db_data/mesh_info.json")
mesh_info2 <- fromJSON(file="data/db_data/mesh_info2.json")
mesh_info3 <- fromJSON(file="data/db_data/mesh_info3.json")

# extract the mappings
mesh_df <- vec_to_df(mesh_info) %>% rename("MeSH"="id")
mesh_df2 <- vec_to_df(mesh_info2) %>% rename("MeSH"="id")
mesh_df3 <- vec_to_df(mesh_info3) %>% rename("MeSH"="id")


mesh_df_full <- do.call(rbind, list(mesh_df, mesh_df2, mesh_df3))

table(mesh_df_full$name == "Add") # 207 of these....
mesh.failed <- filter(mesh_df_full, name=="Add")$MeSH
mesh.failed2 <- mesh.failed[sapply(mesh.failed, grabMeshType) %in% c("C", "D")] # only 19 are actual MeSH Ids - these do not have anything on the webpages

mesh_df_full <- filter(mesh_df_full, name!="Add")
write.table(mesh_df_full, "data/db_data/mesh_info_df.txt", sep="\t", row.names=FALSE)


