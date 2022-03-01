library(tidyverse)
library(data.table)
library(pbapply)
library(pbmcapply)

args <- commandArgs(trailingOnly = TRUE)
counter <- as.integer(args[1])
ngam <- as.integer(args[2])
nsnp <- as.integer(args[3])
cov <- as.numeric(args[4])

#ngam <- 1000 #as.integer(args[2])
#nsnp <- 30000 #as.integer(args[3])
#cov <- 0.01 #as.numeric(args[4])

#ngam <- 50
#nsnp <- 100000
#cov <- 0.511
overlap_denom <- 2

flist <- list.files("~/workzfs-rmccoy22/kweave23/rhapsodi/genModel", "assess_out_ws_", full.names=TRUE, recursive=TRUE)

flist_df <- flist %>% 
  as.data.frame() %>% `colnames<-`("file") %>% 
  separate(file, into=c("base1", "base2", "base3", "base4", "base5", "base6", "base7", "base8", "desc", "rdata"), sep="/", fill="left", remove=FALSE) %>% 
  separate(desc, into=c("ngam", "nsnp", "cov", "se", "avgr", "mse", "mr"), sep = "_", fill="right", remove=FALSE) %>% 
  separate(rdata, into=c("assess", "out", "ws", "window_length", "od", "overlap_denom", "rs", "rsd" ), sep='_', remove=FALSE) %>% 
  mutate_at(vars("ngam", "nsnp", "cov", "se", "avgr", "mse", "mr", "rsd" ), funs(gsub("[gcsemrRdat]", "",.))) %>% 
  mutate(across(c("ngam", "nsnp", "rsd", "window_length"), as.integer)) %>% 
  mutate(across(c("cov", "se", "avgr", "mse", "mr", "overlap_denom"), as.numeric))

roi <- which(flist_df$ngam == ngam & flist_df$nsnp == nsnp & flist_df$cov == cov)
flist_df_subset <- flist_df[roi,c("ngam", "nsnp", "cov", "overlap_denom", "window_length", "se", "avgr", "rsd")]
file_paths_subset <- flist_df[roi,"file"]
message(nrow(flist_df_subset))
if (nrow(flist_df_subset) != counter){
  message("some are missing")
  flist_df_desired <- data.frame(ngam = rep(ngam, counter), 
                                 nsnp = rep(nsnp, counter),
                                 cov = rep(cov, counter), 
                                 overlap_denom = rep(overlap_denom, counter),
                                 window_length = rep(0, counter),
                                 se = rep(0, counter),
                                 avgr = rep(0, counter),
                                 rsd = rep(0, counter))
  iterval <- 1
  windows <- c(250, 500, 750, seq(1000, nsnp, 500))
  #windows2 <- seq(nsnp/overlap_denom, nsnp, 500)
  seqerrs <- c(0.001, 0.005, 0.05)
  avgrs <- c(0.6, 1, 3)
  rsds <- c(42, 357, 1848)
  for (win in windows){
  #for (win in c(windows, windows2)){
    for (seqerr in seqerrs){
      for (avgr in avgrs){
        for (rsd in rsds){
          flist_df_desired[iterval, "window_length"] <- win
          flist_df_desired[iterval, "se"] <- seqerr
          flist_df_desired[iterval, "avgr"] <- avgr
          flist_df_desired[iterval, "rsd"] <-  rsd
          iterval <- iterval + 1
        }
      }
    }
  }
  
  flist_df_desired <- flist_df_desired %>% 
    mutate(across(c("ngam", "nsnp", "rsd", "window_length"), as.integer)) %>% 
    mutate(across(c("cov", "se", "avgr", "overlap_denom"), as.numeric))

  flist_df_subset_unite <- flist_df_subset %>% unite("all_cols") %>% as.data.frame()
  flist_df_desired_unite <- flist_df_desired %>% unite("all_cols") %>% as.data.frame()
  missing <- anti_join(flist_df_desired_unite, flist_df_subset_unite)
  missing
} else {
  message("all are here")
}

load_metrics <- function(file_name){
  load(file_name)
  to_return <- data.frame(phasing_acc = assess_out$phasing$acc, 
                          phasing_ser = assess_out$phasing$ser,
                          phasing_com = assess_out$phasing$com,
                          phasing_lhs = assess_out$phasing$lhs)
  return(to_return)
}

metrics <- do.call(rbind, mclapply(1:length(file_paths_subset),
                                   function(x) load_metrics(file_paths_subset[x]),
                                   mc.cores = getOption("mc.cores", 8)))


flist_df_subset$nrow <- 1:nrow(flist_df_subset)
metrics$nrow <- 1:nrow(flist_df_subset)
to_save <- full_join(flist_df_subset, metrics, by="nrow")

write.csv(to_save, paste0("/scratch/groups/rmccoy22/kweave23/assess_rhapsodi/vary_window/g",ngam, "_s",nsnp,"_c",cov,"_assessed.csv"))
#write.csv(to_save, "/scratch/groups/rmccoy22/kweave23/assess_rhapsodi/vary_window/g50_s100000_c0.511_assessed.csv")
#write.csv(to_save, "/scratch/groups/rmccoy22/kweave23/assess_rhapsodi/vary_window/g2500_s5000_c0.693_assessed.csv")
#write.csv(to_save, "/scratch/groups/rmccoy22/kweave23/assess_rhapsodi/vary_window/g1000_s30000_c0.01_assessed.csv")
