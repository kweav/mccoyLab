library(tidyverse)

nl_df <- read_csv("neale_lab_lsdsc_resource_20220201.csv")
nl_df <- nl_df %>% 
  separate(phenotype, into=c("phenotypeDescA", "phenotypeDescB", "phenotypeDescC", "phenotypeDescD", "normalized"), sep="_", fill="left", remove=FALSE) %>% 
  dplyr::filter(normalized == "irnt" & sex != "both_sexes")
  
description_df <- nl_df %>% 
  dplyr::select(phenotype, description, sex)
wget_df <- nl_df %>% 
  dplyr::select(ldsc_sumstat_wget)

write_csv(description_df, "nl_sumstat_descriptions.csv")
write.table(wget_df, "nl_sumstat_wgets.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)