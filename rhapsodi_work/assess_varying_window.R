library(rhapsodi)

args <- commandArgs(trailingOnly = TRUE)
num_gametes <- args[1]
num_snps <- args[2]
cov <- args[3]
seq_error <- args[4]
avg_recomb <- args[5]
random_seed <- args[6]
threads <- args[7]
window_size <- args[8]
overlap_denom <- args[9]

file_base <- "~/workzfs-rmccoy22/kweave23/rhapsodi/genModel/"
sim_base <- paste0("g",num_gametes, "_s",num_snps,"_c",cov,"_se",seq_error,"_r",avg_recomb,"/runGen_gam_",num_gametes,"_snp_",num_snps,"_cov_",cov,"_seqerr_",seq_error,"_avgr_",avg_recomb,"_rs_",random_seed)

if (seq_error=="0.005" & avg_recomb=="1"){
  sim_out <- paste0("g",num_gametes,"_s",num_snps,"_c",cov,"_se",seq_error,"_r",avg_recomb,"/assess_out_ws_",window_size,"_od_",overlap_denom,"_rs_",random_seed,".Rdata")
  pred_out <- paste0("g",num_gametes,"_s",num_snps,"_c",cov,"_se",seq_error,"_r",avg_recomb,"/rhapsodi_out_ws_",window_size,"_od_",overlap_denom,"_rs_",random_seed,".Rdata")
} else{
  file_base_add <- "changing_model_params/"
  sim_out_dir_base <- paste0("g",num_gametes,"_s",num_snps,"_c",cov,"_se",seq_error,"_r",avg_recomb,"_mse", seq_error, "_mr", avg_recomb)
  sim_out_dir <- paste0(file_base, file_base_add, sim_out_dir_base)
  if (!dir.exists(sim_out_dir)){
    dir.create(sim_out_dir, showWarnings = FALSE, recursive = TRUE)
  }
  sim_out <- paste0(sim_out_dir, "/assess_out_ws_",window_size,"_od_",overlap_denom,"_rs_",random_seed,".Rdata")
  pred_out <- paste0(sim_out_dir, "/rhapsodi_out_ws_",window_size,"_od_",overlap_denom,"_rs_",random_seed,".Rdata")
}

input_file <- paste0(file_base, sim_base, "_gametedf_na_truth_afseqednm.csv")
dt <- read.delim(input_file, sep=",", na.strings = c("NA"))

rhapsodi_out <- rhapsodi::rhapsodi_autorun(NULL, use_dt = TRUE, input_dt = dt,
  mcstop = FALSE, threads=as.integer(threads),
  window_length = as.integer(window_size), overlap_denom = as.numeric(overlap_denom),
  seqError_model = as.numeric(seq_error), avg_recomb_model = as.numeric(avg_recomb),
  sampleName = "simWL", chrom="chrT", verbose=TRUE)

save(rhapsodi_out, file=pred_out)

true_donor_haps_file <- paste0(file_base, sim_base, "_donorHaps_truth_afseqednm.csv")
true_dh <- read.delim(true_donor_haps_file, sep=",", na.strings = c("NA"))

true_gamete_full_file <- paste0(file_base, sim_base, "_gametedf_full_truth_afseqednm.csv")
true_gf <- read.delim(true_gamete_full_file, sep = ",", na.strings = c("NA"))

true_ci_file <- paste0(file_base, sim_base, "_crossoverIndices_truth_ptfseqednm.csv")
true_ci <- read.delim(true_ci_file, sep = ",", na.strings = c("NA"))

assess_out <- list()

assess_phasing_out <- tryCatch(rhapsodi::sim_assess_phasing(true_dh, rhapsodi_out$donor_haps, nrow(true_dh)),
                                error=function(e) {return(2)})
assess_out$phasing <- assess_phasing_out

assess_gam_out <- tryCatch(rhapsodi::sim_assess_gam_imputation(true_gf, rhapsodi_out$unsmoothed_gamete_genotypes, nrow(true_gf), ncol(true_gf[,-1])),
                            error=function(e) {return(2)})
assess_out$gam_imputation <- assess_gam_out

assess_recomb_out <- tryCatch(rhapsodi::sim_assess_recomb(true_ci, rhapsodi_out$recomb_breaks),
                                error=function(e) {return(2)})
assess_out$recomb <- assess_recomb_out

save(assess_out, file=sim_out)
