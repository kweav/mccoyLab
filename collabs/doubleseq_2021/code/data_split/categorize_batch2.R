library(tidyverse)
library(data.table)

random_seed <- 7319
set.seed(random_seed)

embryoSet_file <- 'data_split/embryoID_by_set.csv'
categorized_meta_file <- 'tidied_meta/tidied_meta_batch1_CREATE_kw.csv'
tobe_categorized_meta_file <- 'tidied_meta/tidied_meta_batch2_CREATE_kw.csv'
meta_full_file <- 'tidied_meta/tidied_meta_CREATE_kw_20220308.csv'

embryoSet <- read.csv(embryoSet_file, row.names  = 1)
categorized_meta <- read.csv(categorized_meta_file, row.names=1)
tobe_categorized_meta <- read.csv(tobe_categorized_meta_file)
full_meta <- read.csv(meta_full_file, row.names = 1)

categorizedSet <- merge(embryoSet, categorized_meta, by=0)
study_participant_set <- categorizedSet[,c("set", "Study_Participant_ID")] %>% 
  group_by(Study_Participant_ID) %>% unique


categorize_batch2_pt1 <- merge(study_participant_set, tobe_categorized_meta, by = "Study_Participant_ID", all.y = TRUE) %>% `rownames<-`(.[,"X"]) %>% select(-X)
tobe_still <- categorize_batch2_pt1[which(is.na(categorize_batch2_pt1$set)), ] %>% select(-set) 
categorize_batch2_pt1 <- drop_na(categorize_batch2_pt1, "set")
categorize_batch2_pt1 <- cbind(embryoID = row.names(categorize_batch2_pt1), categorize_batch2_pt1)

###split meta data in a study participant aware manner, limiting a single study participant to only one set
unique_participant_ID <- as.numeric(names(table(tobe_still$Study_Participant_ID)))
counts_participant_ID <- as.numeric(table(tobe_still$Study_Participant_ID))
o_age_by_pid <- unlist(lapply(1:length(unique_participant_ID), function(x) as.numeric(names(table(tobe_still$Oocyte_Age[which(tobe_still$Study_Participant_ID == unique_participant_ID[x])]))))) %>% `names<-`(unique_participant_ID)
s_age_by_pid <- unlist(lapply(1:length(unique_participant_ID), function(x) as.numeric(names(table(tobe_still$Sperm_Age[which(tobe_still$Study_Participant_ID == unique_participant_ID[x])]))))) %>% `names<-`(unique_participant_ID)


### stackoverflow method https://stackoverflow.com/questions/54566428/caret-creating-stratified-data-sets-based-on-several-variables
df <- data.frame("upid"=unique_participant_ID, "embryo_counts"=counts_participant_ID, "oocyte_age" = o_age_by_pid, "sperm_age" = s_age_by_pid) %>% 
  mutate(n = row_number()) %>% 
  dplyr::select(n, everything())
train <- df %>% 
  group_by("embryo_counts", "oocyte_age") %>% 
  sample_frac(0.75)
test <- anti_join(df, train)

train_study_participants <- train$upid
train_embryoID <- rownames(tobe_still)[which(tobe_still$Study_Participant_ID %in% train_study_participants)]
test_study_participants <- test$upid
test_embryoID <- rownames(tobe_still)[which(tobe_still$Study_Participant_ID %in% test_study_participants)]

train_eid_df <- data.frame(embryoID = train_embryoID, set = "train")
test_eid_df <- data.frame(embryoID = test_embryoID, set = "test")
eid_df <- rbind(train_eid_df, test_eid_df) %>% `rownames<-`(.[,"embryoID"])
eid_df <- rbind(eid_df, categorize_batch2_pt1[,c("embryoID", "set")])
#write.csv(eid_df, 'data_split/embryo_bySet_batch2.csv')

full_eid_df <- rbind(embryoSet, select(eid_df, "set"))
#write.csv(full_eid_df, 'data_split/embryo_bySet_full_kw_20211217.csv')
full_meta_withSet <- merge(full_eid_df, full_meta, by=0) %>% `rownames<-`(.[,"Row.names"]) %>% select(-"Row.names")
#write.csv(full_meta_withSet, 'tidied_meta/meta_withSet_kw_20211217.csv')
write.csv(full_meta_withSet, 'tidied_meta/meta_withSet_kw_20220308.csv')