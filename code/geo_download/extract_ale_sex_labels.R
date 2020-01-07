
# Reformat ALE data to extract GSE, GSM, and sex labels.

require('tidyverse')

reform_ale <- function(ale){
 ale %>%
  rename("gsm"=X1, "gpl"=PlatformID, "gse"=ExperimentID, "sex"=Gender) %>% 
  select(-Age, -TissueID, -DiseaseState, -TissueName, -TaxonID, -Molecule) %>%
  mutate(gsm=paste("GSM", gsm, sep=""), 
         gse=paste("GSE", gse, sep=""), 
         gpl=paste("GPL", gpl, sep="")) %>%
  mutate(sex = case_when(
    sex == "M" ~ "male",
    sex == "F" ~ "female"))
}

ale_human <- read_csv("../data/heuristic/ale_out")
ale_human2 <- reform_ale(ale_human) 
ale_human2 %>% write_csv("data/ale_labels/ale_human_reform.csv") 

ale_mouse <- read_csv("../data/heuristic/ale_out_10090")
ale_mouse2 <- reform_ale(ale_mouse) 
ale_mouse2 %>% write_csv("data/ale_labels/ale_mouse_reform.csv") 

ale_rat <- read_csv("../data/heuristic/ale_out_10116")
ale_rat2 <- reform_ale(ale_rat) 
ale_rat2 %>% write_csv("data/ale_labels/ale_rat_reform.csv") 
