library(tidyverse)
library(data.table)

random_seed <-  7319 #5281 #6133 #8888 #9746  #2999 #4263 #7319 #277 #198 #42 #3854
set.seed(random_seed)
file_meta <- "tidied_meta/tidied_meta_batch1_CREATE_kw.csv"

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
write_csv(eid_df, "data_split/embryoID_by_set.csv")

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

meta_count_train_indices <- which(meta$Study_Participant_ID %in% train_study_participants)
trainbar.x <- as.numeric(table(train$embryo_counts))
names(trainbar.x) <- names(table(train$embryo_counts))
barplot(trainbar.x, xlab='Number of Embryos', ylab="Number of Study Participants", cex.names=0.5) + title("Train")

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