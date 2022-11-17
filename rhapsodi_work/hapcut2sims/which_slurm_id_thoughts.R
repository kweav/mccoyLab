library(tidyverse)
library(dplyr)

conditions <- read.delim("~/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/conditions.txt", header=FALSE, col.names = c("ngam", "nsnp", "cov", "seqerr", "avgr", "rsd")) %>%
  mutate(row_id=dplyr::row_number())
slurm_ids <- c(9544995, 9737716, 9761424, 9861373, 9939379)

conditions$which_slurm_file <- 0


new_index_to_row_id <- read.delim("/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/rerun_ids.txt", header = FALSE, col.names = "row_id")
rerun_index_to_row_id <- read.delim("/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/final_rerun_ids.txt", header = FALSE, col.names = "row_id")

for (i in 1:nrow(conditions)){
  which_slurm_file <- 0
  if (file_test("-f", paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[5], "_", which(rerun_index_to_row_id$row_id == i), ".out"))){
    which_slurm_file <- paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[5], "_", which(rerun_index_to_row_id$row_id == i), ".out")
  } else if(file_test("-f", paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[4], "_", which(new_index_to_row_id$row_id == i), ".out"))){
    which_slurm_file <- paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[4], "_", which(new_index_to_row_id$row_id == i), ".out")
  } else if(file_test("-f", paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[1], "_", conditions[i, "row_id"], ".out"))){
    which_slurm_file <- paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[1], "_", conditions[i, "row_id"], ".out")
  } else if (file_test("-f", paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[3], "_", conditions[i, "row_id"], ".out"))){
    which_slurm_file <- paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[3], "_", conditions[i, "row_id"], ".out")
  } else if (file_test("-f", paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[2], "_", conditions[i, "row_id"], ".out"))){
    which_slurm_file <- paste0('/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/slurm-', slurm_ids[2], "_", conditions[i, "row_id"], ".out")
  }
  if (which_slurm_file != 0){
  conditions$which_slurm_file <- which_slurm_file
  } else{
    message("no files")
  }
}

write.csv(conditions, file="/home/kweave23/scr4_rmccoy22/kweave23/wslurmid_conditions.txt")