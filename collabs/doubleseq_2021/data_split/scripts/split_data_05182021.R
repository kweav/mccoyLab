library(DESeq2)
library(SummarizedExperiment)
library(rtracklayer)
library(tidyverse)
library(data.table)
library(caret)
library(splitstackshape)
library(readxl)
library(janitor)
library(ggrepel)
library(ggthemes)
library(ggplot2)
library(edgeR)
library(plotly)
library(survcomp)
library(cowplot)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(VennDiagram)

#`%notin%` <- Negate(`%in%`)

###Files
args <- commandArgs(trailingOnly = TRUE)
where_run <- "local" #args[1]
random_seed <-  7319 #5281 #6133 #8888 #9746  #2999 #4263 #7319 #277 #198 #42 #3854 #as.integer(args[2])
set.seed(random_seed)
if (where_run == "comp") {
  load("/home/kweave23/create/create_summarized_experiment.Rdata")
  file_gencode <- "/home/kweave23/create/resources/gencode.v34.annotation.gtf"
  file_meta <- "/home/kweave23/create/tidied_meta_CREATE_kw.csv"
} else if (where_run == "local"){
  load("/Users/kateweaver/mccoyLab/collabs/doubleseq_2021/create_summarized_experiment.Rdata")
  file_gencode <- "/Users/kateweaver/genomes/hg38_genome/gencode.v34.annotation.gtf"
  file_meta <- "/Users/kateweaver/mccoyLab/collabs/doubleseq_2021/tidied_meta_CREATE_kw.csv"
} else {stop("no known file path to summarized experiment")}

###load meta data
meta <- read.csv(file_meta, row.names=1) %>%
  as.data.frame %>%
  mutate(across(c("AOD", "GC", "Infertility_type", "Previous_pregnancy", "Past_surgical_hist", "Pregnant", "Ongoing_pregnancy", "Final_outcome", "Embryo_grade_at_freezing", "Interpretation", "cDNA_RT_Date", "Library_Prep_Date", "Sequencing_Date", "Study_Participant_ID"), as.factor)) %>%
  mutate(across(c("InfD_SSM_GC", "InfD_Egg_factor", "InfD_MF", "InfD_Uterine_factor", "InfD_TF", "InfD_RPL", "InfD_RIF", "InfD_Unexplained", "PMdH_none", "PMdH_vasculitis", "PMdH_immune", "PMdH_stress_hormones"), as.factor))

###split meta data in a study participant aware manner, limiting a single study participant to only one set
unique_participant_ID <- as.numeric(names(table(meta$Study_Participant_ID)))
counts_participant_ID <- as.numeric(table(meta$Study_Participant_ID))
o_age_by_pid <- c()
s_age_by_pid <- c()
for (upid in unique_participant_ID){
  o_age_by_pid <- c(o_age_by_pid, as.numeric(names(table(meta$Oocyte_Age[which(meta$Study_Participant_ID == upid)]))))
  s_age_by_pid <- c(s_age_by_pid, as.numeric(names(table(meta$Sperm_Age[which(meta$Study_Participant_ID == upid)]))))
}

### stackoverflow method https://stackoverflow.com/questions/54566428/caret-creating-stratified-data-sets-based-on-several-variables
df <- data.frame("upid"=unique_participant_ID, "embryo_counts"=counts_participant_ID, "oocyte_age" = o_age_by_pid, "sperm_age" = s_age_by_pid) %>% 
  mutate(n = row_number()) %>% 
  dplyr::select(n, everything())
train <- df %>% 
  group_by("embryo_counts", "oocyte_age") %>% 
  sample_frac(0.75)
test <- anti_join(df, train)
train_study_participants <- train$upid
train_embryoID <- rownames(meta)[which(meta$Study_Participant_ID %in% train_study_participants)]
test_study_participants <- test$upid
test_embryoID <- rownames(meta)[which(meta$Study_Participant_ID %in% test_study_participants)]

train_eid_df <- data.frame(embryoID = train_embryoID, set = "train")
test_eid_df <- data.frame(embryoID = test_embryoID, set = "test")
eid_df <- rbind(train_eid_df, test_eid_df)
write_csv(eid_df, "embryoID_by_set.csv")

pie_plot <- function(slices, lbls, main_title){
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  pie(slices,labels = lbls, col=rainbow(length(lbls)),
      main=main_title) 
}

violin_plot <- function(data_frame_toplot, ylabel){
  ggplot(data_frame_toplot, aes(x=name, y=val, fill=name)) + scale_x_discrete(limits=c("Train", "Test")) + geom_violin(scale="width", adjust=1, width=0.5) + labs(y=ylabel, x="set") + geom_boxplot(width = 0.2)
}

slices <- c(length(train_embryoID), length(test_embryoID))
lbls <- c("In Train", "In Test")
pie_plot(slices, lbls, "Embryos -- study participant aware split")

slices <- c(length(train_study_participants), length(test_study_participants))
lbls <- c("In Train", "In Test")
pie_plot(slices, lbls, "Study Participants -- study participant aware split")


###

# ### stratified method
# df <- data.frame("upid"=unique_participant_ID, "embryo_counts"=counts_participant_ID, "oocyte_age" = o_age_by_pid, "sperm_age" = s_age_by_pid)
# test_df <- stratified(df, c("oocyte_age", "embryo_counts"), .25)
# test.index <- which(unique_participant_ID %in% test_df$upid)

# train.index <- which(unique_participant_ID %notin% test_df$upid)
# ###

# ### Original Method
# set.seed(42)
# train.index <- caret::createDataPartition(counts_participant_ID, p=0.75, list=FALSE, times=1)
# ###

meta_count_train_indices <- which(meta$Study_Participant_ID %in% train_study_participants)
trainbar.x <- as.numeric(table(train$embryo_counts))
names(trainbar.x) <- names(table(train$embryo_counts))
barplot(trainbar.x, xlab='Number of Embryos', ylab="Number of Study Participants", cex.names=0.5) + title("Train")

# ### Original Method
# test_study_participants <- unique_participant_ID[-train.index]
# ###

# ### stratified method
# test_study_participants <- test_df$upid
# ###

meta_count_test_indices <- which(meta$Study_Participant_ID %in% test_study_participants)
testbar.x <- as.numeric(table(test$embryo_counts))
names(testbar.x) <- names(table(test$embryo_counts))
barplot(testbar.x, xlab='Number of Embryos', ylab="Number of Study Participants", cex.names=0.5) + title("Test")

meta_train <- meta[meta_count_train_indices,]
meta_test <- meta[meta_count_test_indices,]

# Oocyte_Age
train_oa <- data.frame(val=meta_train$Oocyte_Age, name="Train")
test_oa <- data.frame(val=meta_test$Oocyte_Age, name="Test")
oaDat <- rbind(train_oa, test_oa)
violin_plot(oaDat, "Oocyte Age")

# Uterus_Age
train_ua <- data.frame(val=meta_train$Uterus_Age, name='Train')
test_ua <- data.frame(val=meta_test$Uterus_Age, name='Test')
uaDat <- rbind(train_ua, test_ua)
violin_plot(uaDat, "Uterus Age")

# Sperm_Age
train_sa <- data.frame(val=meta_train$Sperm_Age, name='Train')
test_sa <- data.frame(val=meta_test$Sperm_Age, name='Test')
saDat <- rbind(train_sa, test_sa)
violin_plot(saDat, "Sperm Age")

# BMI
train_bmi <- data.frame(val=meta_train$BMI, name='Train')
test_bmi <- data.frame(val=meta_test$BMI, name='Test')
bmiDat <- rbind(train_bmi, test_bmi)
violin_plot(bmiDat, "BMI")

# Embryo_grade_at_freezing
train_egaf <- data.frame(val=meta_train$Embryo_grade_at_freezing, name='Train')
test_egaf <- data.frame(val=meta_test$Embryo_grade_at_freezing, name='Test')
egafDat <- rbind(train_egaf, test_egaf)

train_egaf_x <- as.numeric(table(train_egaf$val))
names(train_egaf_x) <- names(table(train_egaf$val))
barplot(train_egaf_x, cex.names=0.35, xlab="Embryo Grade at Freezing", ylab="Number of embryos") + title("Train")

test_egaf_x <- as.numeric(table(test_egaf$val))
names(test_egaf_x) <- names(table(test_egaf$val))
barplot(test_egaf_x, cex.names=0.35, xlab="Embryo Grade at Freezing", ylab="Number of embryos") + title("Test")

# lining_thickness_mm
train_lt <- data.frame(val=meta_train$lining_thickness_mm, name='Train')
test_lt <- data.frame(val=meta_test$lining_thickness_mm, name='Test')
ltDat <- rbind(train_lt, test_lt)
violin_plot(ltDat, "Lining thickness (mm)")

# Pregnant
train_p <- data.frame(val=meta_train$Pregnant, name='Train')
test_p <- data.frame(val=meta_test$Pregnant, name='Test')
pDat <- rbind(train_p, test_p)

train_p_x <- as.numeric(table(train_p$val))
names(train_p_x) <- c("no", "yes")
barplot(train_p_x, cex.names=0.5, xlab="Pregnant", ylab="Number of embryos") + title("Train")

test_p_x <- as.numeric(table(test_p$val))
names(test_p_x) <- c("no", "yes")
barplot(test_p_x, cex.names=0.5, xlab="Pregnant", ylab="Number of embryos") + title("Test")

# cDNA_RT_Date
train_cdna <- data.frame(val=meta_train$cDNA_RT_Date, name='Train')
test_cdna <- data.frame(val=meta_test$cDNA_RT_Date, name='Test')
cdnaDat <- rbind(train_cdna, test_cdna)

train_cdna_x <- as.numeric(table(train_cdna$val))
names(train_cdna_x) <- names(table(train_cdna$val))
barplot(train_cdna_x, cex.names=0.35, xlab="cDNA RT Batch", ylab="Number of embryos") + title("Train")

test_cdna_x <- as.numeric(table(test_cdna$val))
names(test_cdna_x) <- names(table(test_cdna$val))
barplot(test_cdna_x, cex.names=0.35, xlab="cDNA RT Batch", ylab="Number of embryos") + title("Test")

# Library_Prep_Date
train_lpd <- data.frame(val=meta_train$Library_Prep_Date, name='Train')
test_lpd <- data.frame(val=meta_test$Library_Prep_Date, name='Test')
lpdDat <- rbind(train_lpd, test_lpd)

train_lpd_x <- as.numeric(table(train_lpd$val))
names(train_lpd_x) <- names(table(train_lpd$val))
barplot(train_lpd_x, cex.names=0.5, xlab="Library Prep Batch", ylab="Number of embryos") + title("Train")

test_lpd_x <- as.numeric(table(test_lpd$val))
names(test_lpd_x) <- names(table(test_lpd$val))
barplot(test_lpd_x, cex.names=0.5, xlab="Library Prep Batch", ylab="Number of embryos") + title("Test")


# Sequencing_Date
train_sd <- data.frame(val=meta_train$Sequencing_Date, name='Train')
test_sd <- data.frame(val=meta_test$Sequencing_Date, name='Test')
sdDat <- rbind(train_sd, test_sd)

train_sd_x <- as.numeric(table(train_sd$val))
names(train_sd_x) <- names(table(train_sd$val))
barplot(train_sd_x, cex.names=0.5, xlab="Sequencing Batch", ylab="Number of embryos") + title("Train")

test_sd_x <- as.numeric(table(test_sd$val))
names(test_sd_x) <- names(table(test_sd$val))
barplot(test_sd_x, cex.names=0.5, xlab="Sequencing Batch", ylab="Number of embryos") + title("Test")

###split meta data in a random manner, allowing a single study participant to be in more than one set
meta %>%
  mutate(n = row_number()) %>%
  dplyr::select(n, everything())
meta$Embryo_ID <- rownames(meta)
train2 <- meta %>%
  group_by("Oocyte_Age", "Sperm_age") %>%
  sample_frac(0.75)
test2 <- anti_join(meta, train2)
train_study_participants2 <- unique(train2$Study_Participant_ID)
test_study_participants2 <- unique(test2$Study_Participant_ID)
sp_df <- merge(data.frame(study_participant = sort(train_study_participants2), cohort=rep("train", length(train_study_participants2))),
               data.frame(study_participant = sort(test_study_participants2), cohort=rep("test", length(test_study_participants2))), by=1)
shared_sp <- sp_df$study_participant
num_shared_sp <- length(shared_sp)
num_embryos_from_shared_sp <- length(which(meta$Study_Participant_ID %in% shared_sp))

slices <- c(length(train_embryoID2), length(test_embryoID2))
lbls <- c("In Train", "In Test")
pie_plot(slices, lbls, "Embryos -- random split")

slices <- c(length(train_study_participants2) - num_shared_sp, length(test_study_participants2) - num_shared_sp, num_shared_sp)
lbls <- c("Unique to Train", "Unique to Test", "Shared")
pie_plot(slices, lbls,"Study Participants")

slices <- c(length(train_study_participants2) - num_shared_sp, num_shared_sp)
lbls <- c("Unique to Train", "Shared with Test")
pie_plot(slices, lbls,"Train Study Participants")

slices <- c(length(test_study_participants2) - num_shared_sp, num_shared_sp)
lbls <- c("Unique to Test", "Shared with Train")
pie_plot(slices, lbls,"Test Study Participants")

slices <- c(nrow(meta) - num_embryos_from_shared_sp, num_embryos_from_shared_sp)
lbls <- c("From Unshared Study Participants", "From Shared Study Participants")
pie_plot(slices, lbls,"Embryos")

train_embryoID2 <- train2$Embryo_ID
test_embryoID2 <- test2$Embryo_ID

#Oocyte_Age
train_oa2 <- data.frame(val=train2$Oocyte_Age, name="Train")
test_oa2 <- data.frame(val=test2$Oocyte_Age, name="Test")
oaDat2 <- rbind(train_oa2, test_oa2)
violin_plot(oaDat2, "Oocyte Age")

#Uterus_Age
train_ua2 <- data.frame(val=train2$Uterus_Age, name='Train')
test_ua2 <- data.frame(val=test2$Uterus_Age, name='Test')
uaDat2 <- rbind(train_ua2, test_ua2)
violin_plot(uaDat2, "Uterus Age")

#Sperm_Age
train_sa2 <- data.frame(val=train2$Sperm_Age, name='Train')
test_sa2 <- data.frame(val=test2$Sperm_Age, name='Test')
saDat2 <- rbind(train_sa2, test_sa2)
violin_plot(saDat2, "Sperm Age")

# BMI
train_bmi2 <- data.frame(val=train2$BMI, name='Train')
test_bmi2 <- data.frame(val=test2$BMI, name='Test')
bmiDat2 <- rbind(train_bmi2, test_bmi2)
violin_plot(bmiDat2, "BMI")

# Embryo_grade_at_freezing
train_egaf2 <- data.frame(val=train2$Embryo_grade_at_freezing, name='Train')
test_egaf2 <- data.frame(val=test2$Embryo_grade_at_freezing, name='Test')
egafDat2 <- rbind(train_egaf2, test_egaf2)

train_egaf_x2 <- as.numeric(table(train_egaf2$val))
names(train_egaf_x2) <- names(table(train_egaf2$val))
barplot(train_egaf_x2, cex.names=0.35, xlab="Embryo Grade at Freezing", ylab="Number of embryos") + title("Train")

test_egaf_x2 <- as.numeric(table(test_egaf2$val))
names(test_egaf_x2) <- names(table(test_egaf2$val))
barplot(test_egaf_x2, cex.names=0.35, xlab="Embryo Grade at Freezing", ylab="Number of embryos") + title("Test")

# lining_thickness_mm
train_lt2 <- data.frame(val=train2$lining_thickness_mm, name='Train')
test_lt2 <- data.frame(val=test2$lining_thickness_mm, name='Test')
ltDat2 <- rbind(train_lt2, test_lt2)
violin_plot(ltDat2, "Lining thickness (mm)")

#Pregnant
train_p2 <- data.frame(val=train2$Pregnant, name='Train')
test_p2 <- data.frame(val=test2$Pregnant, name='Test')
pDat2 <- rbind(train_p2, test_p2)

train_p_x2 <- as.numeric(table(train_p2$val))
names(train_p_x2) <- c("no", "yes")
barplot(train_p_x2, cex.names=0.5, xlab="Pregnant", ylab="Number of embryos") + title("Train")

test_p_x2 <- as.numeric(table(test_p2$val))
names(test_p_x2) <- c("no", "yes")
barplot(test_p_x2, cex.names=0.5, xlab="Pregnant", ylab="Number of embryos") + title("Test")


#cDNA_RT_Date
train_cdna2 <- data.frame(val=train2$cDNA_RT_Date, name='Train')
test_cdna2 <- data.frame(val=test2$cDNA_RT_Date, name='Test')
cdnaDat2 <- rbind(train_cdna2, test_cdna2)

train_cdna_x2 <- as.numeric(table(train_cdna2$val))
names(train_cdna_x2) <- names(table(train_cdna2$val))
barplot(train_cdna_x2, cex.names=0.35, xlab="cDNA RT Batch", ylab="Number of embryos") + title("Train")

test_cdna_x2 <- as.numeric(table(test_cdna2$val))
names(test_cdna_x2) <- names(table(test_cdna2$val))
barplot(test_cdna_x2, cex.names=0.35, xlab="cDNA RT Batch", ylab="Number of embryos") + title("Test")

#Library_Prep_Date
train_lpd2 <- data.frame(val=train2$Library_Prep_Date, name='Train')
test_lpd2 <- data.frame(val=test2$Library_Prep_Date, name='Test')
lpdDat2 <- rbind(train_lpd2, test_lpd2)

train_lpd_x2 <- as.numeric(table(train_lpd2$val))
names(train_lpd_x2) <- names(table(train_lpd2$val))
barplot(train_lpd_x2, cex.names=0.5, xlab="Library Prep Batch", ylab="Number of embryos") + title("Train")

test_lpd_x2 <- as.numeric(table(test_lpd2$val))
names(test_lpd_x2) <- names(table(test_lpd2$val))
barplot(test_lpd_x2, cex.names=0.5, xlab="Library Prep Batch", ylab="Number of embryos") + title("Test")


#Sequencing_Date
train_sd2 <- data.frame(val=train2$Sequencing_Date, name='Train')
test_sd2 <- data.frame(val=test2$Sequencing_Date, name='Test')
sdDat2 <- rbind(train_sd2, test_sd2)

train_sd_x2 <- as.numeric(table(train_sd2$val))
names(train_sd_x2) <- names(table(train_sd2$val))
barplot(train_sd_x2, cex.names=0.5, xlab="Sequencing Batch", ylab="Number of embryos") + title("Train")

test_sd_x2 <- as.numeric(table(test_sd2$val))
names(test_sd_x2) <- names(table(test_sd2$val))
barplot(test_sd_x2, cex.names=0.5, xlab="Sequencing Batch", ylab="Number of embryos") + title("Test")
