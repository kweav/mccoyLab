library(tidyverse)
library(dplyr)

conditions <- read.delim("~/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/conditions.txt", header=FALSE, col.names = c("ngam", "nsnp", "cov", "seqerr", "avgr", "rsd")) %>%
  mutate(row_id=dplyr::row_number())

#final batch job is 9939379 and is using the `final_rerun_ids.txt` file and fixed the writefmf error and 4 memory not allocated

#current batch job is 9861373 and is just running the IDs from memory allocation or no known reason in `rerun_ids.txt`

#SLURM numbers I want
#9761424 -- latest run, covers those whose directories were not made
#9756074 -- middle run with mkdir -p, but no new IDs from this
#9737716 -- earlyish run
#9544995 -- earliest run, should be most of the <150 gamete ones
slurm_ids <- c(9544995, 9737716, 9761424, 9861373, 9939379)

conditions$hapcut_success <- FALSE
conditions$no_het_snps <- FALSE
conditions$myn_arg <- FALSE
conditions$timeout <- FALSE
conditions$memory_not_allocated <- FALSE
conditions$killed <- FALSE
conditions$hapcut_mad <- FALSE

new_index_to_row_id <- read.delim("/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/rerun_ids.txt", header = FALSE, col.names = "row_id")
rerun_index_to_row_id <- read.delim("/home/kweave23/scr4_rmccoy22/kweave23/rhapsodi_hapcut2_sims/final_rerun_ids.txt", header = FALSE, col.names = "row_id")

for (i in 1:nrow(conditions)){
  if (conditions[i, "seqerr"] == 0.005 & conditions[i, "avgr"] == 1){
    pathoi <- paste0("/home/kweave23/data_rmccoy22/sperm_seq_rhapsodi/genModel_output/mse_0.005_mar_1/g",conditions[i, "ngam"],"_s", as.character(as.integer(conditions[i,"nsnp"])), "_c", conditions[i, "cov"], "_se", conditions[i, "seqerr"], "_r", conditions[i, "avgr"])
  } else {
    pathoi <- paste0("/home/kweave23/data_rmccoy22/sperm_seq_rhapsodi/genModel_output/changing_model_params/g",conditions[i, "ngam"],"_s", as.character(as.integer(conditions[i,"nsnp"])), "_c", conditions[i, "cov"], "_se", conditions[i, "seqerr"], "_r", conditions[i, "avgr"], "_mse", conditions[i, "seqerr"], "_mr", conditions[i, "avgr"])
  }
  hapcut2_outfile_oi <- paste0("rs", conditions[i, "rsd"], "_hapcut2.phased.VCF")
  hapcut_file_test <- file_test("-f",paste0(pathoi, "/", list.files(path = pathoi, hapcut2_outfile_oi)))
  conditions$hapcut_success[i]<- hapcut_file_test
  
  if (!hapcut_file_test){
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
    hetsnpgrep <- grep("no hetSNPs so conversion is exiting", readLines(which_slurm_file), value=TRUE)
    if (!identical(hetsnpgrep, character(0))){
      conditions$no_het_snps[i] <- TRUE
    } else { 
      myngrep <- grep("No fmf conversion possible due to myn argument", readLines(which_slurm_file), value = TRUE)
      if (!identical(myngrep, character(0))){
        conditions$myn_arg[i] <- TRUE
      } else {
        timeoutgrep <- grep("TIME LIMIT", readLines(which_slurm_file), value = TRUE)
          if (!identical(timeoutgrep, character(0))){
            conditions$timeout[i] <- TRUE
          } else {
            memory_alloc_grep <- grep("Cannot allocate memory", readLines(which_slurm_file), value = TRUE)
            if (!identical(memory_alloc_grep, character(0))){
              conditions$memory_not_allocated[i] <- TRUE
            } else{
              hapcut_mad_grep <- grep("No output VCF", readLines(which_slurm_file), value = TRUE)
              if (!identical(hapcut_mad_grep, character(0))){
                conditions$hapcut_mad[i] <- TRUE
              } else {
                killed_grep <- grep("Killed", readLines(which_slurm_file), value = TRUE)
                if (!identical(killed_grep, character(0))){
                  conditions$killed[i] <- TRUE
                }
              }
            }
          }
        }
    }
  }
}

subset_succeeded <- conditions[which(conditions$hapcut_success == TRUE),]
conditions <- conditions[which(conditions$hapcut_success == FALSE),]

subset_timedout <- conditions[which(conditions$timeout == TRUE),]
conditions <- conditions[which(conditions$timeout == FALSE),]

subset_no_hetsnps <- conditions[which(conditions$no_het_snps == TRUE),]
conditions <- conditions[which(conditions$no_het_snps == FALSE),]

subset_myn_arg <- conditions[which(conditions$myn_arg == TRUE),]
conditions <- conditions[which(conditions$myn_arg == FALSE),]

subset_killed <- conditions[which(conditions$killed == TRUE),]
conditions <- conditions[which(conditions$killed == FALSE),]

subset_memory_not_allocated <- conditions[which(conditions$memory_not_allocated == TRUE),]
conditions <- conditions[which(conditions$memory_not_allocated == FALSE),]

subset_hapcut_mad <- conditions[which(conditions$hapcut_mad == TRUE),]
conditions <- conditions[which(conditions$hapcut_mad == FALSE),]


#14 is no VCF but FMF generated 3 30000 0.01  0.001    3  42    277; it's got that thing where the info for the line is instead multiline... Let's play with this file locally.
#> hapcutmad_rowids
#[1]  277  279  286  288  542  545 1216 1221 1222 1224 1225 1230 1231 1233
#[15] 1237 1239 1240 1242 1466 1467 1475 1476 1478 1481 1484 1702 1711 1720
#[29] 1722 2190 2193 2199 2202 2208

#> conditions[279,]
#     ngam  nsnp  cov seqerr avgr  rsd
#279    3 30000 0.01  0.001    3 1848

#> conditions[286,]
#     ngam  nsnp  cov seqerr avgr rsd
#286    3 30000 0.01  0.005    3  42

#> conditions[288,]
#   ngam  nsnp  cov seqerr avgr  rsd
#288    3 30000 0.01  0.005    3 1848

#> conditions[542,]
#   ngam   nsnp cov seqerr avgr rsd
#542    3 100000 0.1  0.001  0.6 357

#> > conditions[2193,]
#   ngam nsnp   cov seqerr avgr  rsd
#2193  150 5000 0.001  0.001    1 1848
