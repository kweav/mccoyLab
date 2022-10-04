library(vcfR)
library(rhapsodi)
library(tidyverse)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

#vcf_file <- vcfR::read.vcfR("~/mccoyLab/rhapsodi_work/hapcut2sims/hapcut2_out_g1000_s30000_cov_0.001_seqerr_0.005_avgr_1_rs42.phased.VCF")
vcf_file <- vcfR::read.vcfR(args[1])
genotypes <- vcfR::extract.gt(vcf_file) %>% as.data.frame() %>%
  mutate(index=dplyr::row_number()) %>% #create index column
  tibble::rownames_to_column("chr_pos") %>% #rownames to column
  tidyr::separate("chr_pos", c("chr", "pos"), sep='_', remove=TRUE) %>% #separate into chr and pos column
  na_if("0/1") %>%
  tidyr::separate(SAMPLENAME, c("h1", "h2"), sep="\\|", remove=TRUE) %>%
  dplyr::select(c("index", "pos", "h1", "h2")) %>% #select index, pos, h1, h2
  mutate(h1 = as.integer(h1), h2 = as.integer(h2)) #make integers

#true_donor_haps_file <- "~/mccoyLab_withOthers/transmission-distortion/test_data_rhapsodi_gen/test_data/g1000_s30000_c0.001_se0.005_r1/runGen_gam_1000_snp_30000_cov_0.001_seqerr_0.005_avgr_1_rs_42_donorHaps_truth_afseqednm.csv"
true_donor_haps_file <- args[2]
true_dh <- read.delim(true_donor_haps_file, sep=",", na.strings = c("NA"))

if (nrow(true_dh) == nrow(genotypes)){
  assess_out <- rhapsodi::sim_assess_phasing(true_dh, genotypes, nrow(true_dh))
  save(assess_out, file = paste0(args[3], "_hapcut2_assess_out.Rdata"))
  save(genotypes, file=paste0(args[3], "_hapcut2_phased_out.Rdata"))
} else{
  message("want to figure out which positions are missing")
}