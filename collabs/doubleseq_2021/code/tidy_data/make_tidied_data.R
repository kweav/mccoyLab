library(tidyverse)
library(data.table)

meta_raw1 <- "raw_meta/DoubleSeq_copied_20210408_kw.csv"
meta_raw2 <- "raw_meta/DoubleSeqBatch2_20211206.csv"

tidy_raw_metadata <- function(meta_file){ 
  comb_col_names <- c("Study_Participant_ID", "Embryo_ID", "AOD", "GC", "Oocyte_Age",
                      "Uterus_Age", "Sperm_Age", "BMI", "Infertility_type",
                      "Infertility_diagnosis", "Previous_pregnancy","Past_med_hist",
                      "Past_surgical_hist", "FSH_D3", "AMH", "AFC", "LH", "E2_Day_2_3",
                      "Type_of_trigger", "Total_Gn", "E2_on_trigger", "Num_eggs",
                      "Num_MII", "Num_2PN", "Num_Blasts", "TMC", "Time_to_OPU_hr",
                      "Time_to_stripping_hr", "Time_to_ICSI_hr", "Time_to_Biopsy_hr",
                      "Embryo_grade_at_freezing", "Interpretation", "Natural_cycle",
                      "Medication", "lining_thickness_mm", "Prednisone", "Fragmin",
                      "EG", "HCG", "IL", "GCSF", "Metformin","PRP","Sildenafil",
                      "Transferring_physician", "Transfer_catheter_used", "Pregnant", "Ongoing_pregnancy",
                      "Final_outcome", "cDNA_RT_Date", "Library_Prep_Date", "Sequencing_Date")
  
  meta <- read.csv(meta_file, header = FALSE, skip=1, col.names = comb_col_names, row.names = 2) %>% 
    as.data.frame %>% 
    na_if(".") %>% 
    na_if("TBD") %>% 
    mutate(AMH = as.numeric(as.character(recode(AMH, ">127"="128")))) %>%
    mutate(LH = as.numeric(as.character(recode(as.character(LH), ">9"="10")))) %>%
    mutate(E2_Day_2_3 = as.numeric(as.character(recode(E2_Day_2_3, "<18"="17")))) %>%
    mutate(Embryo_grade_at_freezing = recode(Embryo_grade_at_freezing, "1AA"="Good", "1AB"="Good", "1BA"="Fair", "1BB"="Fair", "1CB"="Poor", "2AA"="Good", "2AB"="Good", "2BA"="Fair", "2BB"="Fair", "2CC"="Poor", "3AA"="Good", "3BB"="Fair")) %>%
    mutate(Interpretation = recode(Interpretation, "EUPLOID"=0, "MOSAIC"=1)) %>%
    mutate(Prednisone = as.factor(recode(as.character(Prednisone), "B"="1"))) %>%
    mutate(Transferring_physician = as.factor(recode(as.character(Transferring_physician), "PS" = "1", "AB" = "2", "KG" = "3", "CL"="4", "Fellow"="5", "AK"="6"))) %>%
    mutate(Transfer_catheter_used = as.factor(recode(as.character(Transfer_catheter_used), "TomCat"="1", "Tomcat"="1", "TOMCAT"="1", "Wallace"="2", "Cook"="3"))) %>%
    mutate(cDNA_RT_Date = recode(cDNA_RT_Date,  "19-Mar-20"="1", "13-Apr-20"="2", "17-Jul-20"="3", "12-Aug-20"="4", "18-Aug-20"="5", "25-Aug-20"="6", "2-Sep-20"="7", "8-Sep-20"="8", "23-Sep-20"="9", "16-Oct-20"="10", "27-Oct-20"="11", "5-Nov-20"="12", "9-Dec-20"="13", "14-Jan-21"="14", "22-Feb-21"="15", "26-Mar-21"="16", "5-May-21"="17", "27-May-21"="18", "22-Jun-21"="19", "20-Jul-21"="20")) %>%
    mutate(Library_Prep_Date = recode(Library_Prep_Date, "11-Dec-20"="0", "11-Feb-21"="1", "12-Aug-21"="2")) %>%
    mutate(Sequencing_Date = recode(Sequencing_Date, "14-Dec-20"= "0", "5-Mar-21"="1", "28-Sep-21"="2")) %>%
    mutate(across(c("AOD", "GC", "Infertility_type", "Previous_pregnancy", "Past_surgical_hist", "Pregnant", "Ongoing_pregnancy", "Final_outcome", "Embryo_grade_at_freezing", "Interpretation", "cDNA_RT_Date", "Library_Prep_Date", "Sequencing_Date"), as.factor))
  
  
  base_vals <- c("0", "1", "2", "3", "4", "5", "6", "7")
  new_col <- c()
  base_vals2 <- c("0", "1", "2", "3")
  new_col2 <- c()
  for (i in c(1:nrow(meta))){
    new_vector <- rep(0, 8)
    obs_val <- as.character(meta$Infertility_diagnosis[i])
    for (valOI in base_vals){
      index_val <- as.integer(valOI) + 1
      if (valOI %in% unlist(strsplit(obs_val, ","))){
        new_vector[index_val] <- 1
      }
    }
    new_col <- c(new_col, paste(new_vector, collapse="_"))
    new_vector2 <- rep(0, 4)
    obs_val2 <- as.character(meta$Past_med_hist[i])
    for (valOI2 in base_vals2){
      index_val2 <- as.integer(valOI2) + 1
      if (valOI2 %in% unlist(strsplit(obs_val2, ","))){
        new_vector2[index_val2] <- 1
      }
    }
    new_col2 <- c(new_col2, paste(new_vector2, collapse="_"))
  }
  meta$NInfD <- new_col
  meta$NPMdH <- new_col2
  
  meta <- separate(meta, col="NInfD", into=c("InfD_SSM_GC", "InfD_Egg_factor", "InfD_MF", "InfD_Uterine_factor", "InfD_TF", "InfD_RPL", "InfD_RIF", "InfD_Unexplained"), sep='_') %>% 
    separate(col="NPMdH", into=c("PMdH_none", "PMdH_vasculitis", "PMdH_immune", "PMdH_stress_hormones"), sep="_") %>% 
    mutate(across(c("InfD_SSM_GC", "InfD_Egg_factor", "InfD_MF", "InfD_Uterine_factor", "InfD_TF", "InfD_RPL", "InfD_RIF", "InfD_Unexplained", "PMdH_none", "PMdH_vasculitis", "PMdH_immune", "PMdH_stress_hormones"), as.factor)) %>% 
    within(rm("Infertility_diagnosis", "Past_med_hist"))
  
  return(meta)

}

meta_tidied1 <- tidy_raw_metadata(meta_raw1)
write.csv(meta_tidied1, "tidied_meta/tidied_meta_batch1_CREATE_kw.csv")
meta_tidied2 <- tidy_raw_metadata(meta_raw2)
write.csv(meta_tidied2, "tidied_meta/tidied_meta_batch2_CREATE_kw.csv")

meta_full <- rbind(meta_tidied1,  meta_tidied2)
#write.csv(meta_full, "tidied_meta/tidied_meta_CREATE_kw_20211217.csv")
write.csv(meta_full, "tidied_meta/tidied_meta_CREATE_kw_20220308.csv")

