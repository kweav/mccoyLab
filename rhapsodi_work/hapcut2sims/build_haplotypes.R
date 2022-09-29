library(vcfR)
library(rhapsodi)
library(tidyverse)
library(tidyr)

vcf_file <- vcfR::read.vcfR("~/mccoyLab/rhapsodi_work/hapcut2sims/hapcut2_out_g1000_s30000_cov_0.001_seqerr_0.005_avgr_1_rs42.phased.VCF")
genotypes <- vcfR::extract.gt(vcf_file) %>% as.data.frame() %>% 
  na_if("0/1") %>%
  tidyr::separate(SAMPLENAME, c("hap1", "hap2"), sep="\\|", remove=FALSE)