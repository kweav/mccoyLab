library(limma)
library(SummarizedExperiment)
library(rtracklayer)
library(tidyverse)
library(data.table)
library(readxl)
library(janitor)
library(ggrepel)
library(ggthemes)
library(edgeR)
library(plotly)
library(survcomp)
library(cowplot)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)

#importing SummarizedExperiment Data and transforming to raw counts data
#load("/home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/create_outcome_pred/first_de_analysis/create_summarized_experiment.Rdata")
load("~/mccoyLab/collabs/doubleseq_2021/create_summarized_experiment.Rdata")
#load("/home/kweave23/create/create_summarized_experiment.Rdata")
counts <- as.matrix(assays(se)$counts)
for (col in 1:ncol(counts)){
  colnames(counts)[col]<-strsplit(sub("Aligned.sortedByCoord.out.bam", "", colnames(counts)[col]), "_")[[1]][1]
}

#subset to just training data
#embryoID_by_set <- read.csv("/home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/create_outcome_pred/training_test_split_study_participant_aware/embryoID_by_set.csv")
embryoID_by_set <- read.csv("~/mccoyLab/collabs/doubleseq_2021/data_split/embryoID_by_set.csv")
#embryoID_by_set <- read_csv("/home/kweave23/create/embryoID_by_set.csv")
train_embryoIDs = embryoID_by_set$embryoID[embryoID_by_set$set == "train"]
train_cols <- which(colnames(counts) %in% train_embryoIDs)

counts <- counts[,train_cols]

hist(colSums(counts), breaks=50) #library size is right-skewed, normal distribution, but fairly variable
message(max(colSums(counts))/min(colSums(counts))) #nearly 12-fold difference between largest and smallest library size for training samples before filtering. Should probably use voom before limma-trend

#subset to genes with matchind gencode ensembl gene IDs on chr1-22
#file_gencode "home-3/kweave23@jhu.edu//workzfs-rmccoy22/resources/reference/gencode.v34.annotation.gtf"
file_gencode <- "~/genomes/hg38_genome/gencode.v34.annotation.gtf"
#file_gencode <- "/home/kweave23/create/resources/gencode.v34.annotation.gtf"
gtf <- rtracklayer::import(file_gencode) %>% 
  as.data.frame() %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::select(gene_id, seqnames, width) %>% 
  dplyr::rename(ensembl_gene_id = gene_id) %>% 
  dplyr::rename(chromosome_name = seqnames) %>% 
  dplyr::rename(length = width)

gene_table <- gtf[match(rownames(counts), gtf$ensembl_gene_id),]
gene_table <- gene_table[gene_table$chromosome_name %in% paste0("chr", 1:22),]
counts <- counts[gene_table$ensembl_gene_id,]

#Importing metadata to differentiate between pregnant and not pregnant
#meta <- read.csv("/home/kweave23/create/tidied_meta_CREATE_kw.csv", row.names = 1) %>% 
#meta <- read.csv("/home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/create_outcome_pred/tidied_meta/tidied_meta_CREATE_kw.csv", row.names = 1) %>% 
meta <- read.csv("~/mccoyLab/collabs/doubleseq_2021/tidied_meta_CREATE_kw.csv", row.names = 1) %>% 
  as.data.frame() %>% 
  mutate(across(c("AOD", "GC", "Infertility_type", "Previous_pregnancy", "Past_surgical_hist", "Pregnant", "Ongoing_pregnancy", "Final_outcome", "Embryo_grade_at_freezing", "Interpretation", "cDNA_RT_Date", "Library_Prep_Date", "Sequencing_Date", "Study_Participant_ID"), as.factor)) %>%
  mutate(across(c("InfD_SSM_GC", "InfD_Egg_factor", "InfD_MF", "InfD_Uterine_factor", "InfD_TF", "InfD_RPL", "InfD_RIF", "InfD_Unexplained", "PMdH_none", "PMdH_vasculitis", "PMdH_immune", "PMdH_stress_hormones"), as.factor))

#subsetting meta data to just the training rows            
remove <- setdiff(rownames(meta), colnames(counts))
if (length(which(rownames(meta) %in%  remove)) > 0){
  meta <- meta[-which(rownames(meta) %in%  remove),]
}

#need to sort so name order is the same for setting up design for DESeq experiment
counts_order <- order(colnames(counts))
meta_order <- order(rownames(meta))

countdata <- counts[, counts_order]
coldata <- meta[meta_order, ]

hist(colSums(countdata), breaks=50) #library size is right-skewed, normal distribution, but fairly variable
message(max(colSums(countdata))/min(colSums(countdata))) #nearly 16-fold difference between largest and smallest library size after filtering. Should definitely use voom before limma-trend

#Filtering out counts 0 and 1
keep <- rowSums(countdata) > 1
countdata <- countdata[keep, ]

#set up design matrix for edgeR, limma, & voom
#Response variable / outcome
Pregnant <- factor(meta$Pregnant)
# main covariates to consider
Batch <- factor(meta$Sequencing_Date)
# participant age
Age_O <- as.numeric(meta$Oocyte_Age)
#Age_U <- factor(meta$Uterus_Age)
#lining thickness
lthick <- as.numeric(meta$lining_thickness_mm)
#embryo transfer add_ons
pred_ao <- factor(meta$Prednisone)
frag_ao <- factor(meta$Fragmin)
eg_ao <- factor(meta$EG)
hcg_ao <- factor(meta$HCG)
il_ao <- factor(meta$IL)
gcsf_ao <- factor(meta$GCSF)
met_ao  <- factor(meta$Metformin)
#prp_ao <- factor(meta$PRP) (all 0's)
#sil_ao <- factor(meta$Sildenafil) (all 0's)
#embryo grade at freezing
egaf <- factor(meta$Embryo_grade_at_freezing)
#day of development on biopsy / hours from trigger to OPU
dodob <- as.numeric(meta$Time_to_OPU_hr)
design0 <- model.matrix( ~ Pregnant + Batch)
design1 <- model.matrix( ~ Pregnant + Batch + Age_O)
design2 <- model.matrix( ~ Pregnant + Batch + Age_O + lthick)
design3 <- model.matrix( ~ Pregnant + Batch + Age_O + lthick + egaf)
design4 <- model.matrix( ~ Pregnant + Batch + Age_O + lthick + egaf + dodob)
design5 <- model.matrix( ~ Pregnant + Batch + Age_O + lthick + egaf + dodob
                         + pred_ao + frag_ao + eg_ao  + hcg_ao + il_ao + gcsf_ao + met_ao)

dge <- edgeR::DGEList(counts = countdata, samples = coldata)

##Filtering out low count genes
#keep <- edgeR::filterByExpr(dge, design)
#dge <- dge[keep, , keep.lib.sizes=FALSE]

##TMM normalization
#dge <- edgeR::calcNormFactors(dge)

#log counts per million and use of prior.count to damp down the variance of logs of low counts
logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=3)

#requires statmod package
corfit5 <- duplicateCorrelation(logCPM, design=design5, ndups=1, block=meta$Study_Participant_ID) # A slow computation
message(corfit5$consensus.correlation)
g5 <- boxplot(tanh(corfit5$atanh.correlations))

corfit4 <- duplicateCorrelation(logCPM, design=design4, ndups=1, block=meta$Study_Participant_ID) # A slow computation
message(corfit4$consensus.correlation)
g4 <- boxplot(tanh(corfit4$atanh.correlations))

corfit3 <- duplicateCorrelation(logCPM, design=design3, ndups=1, block=meta$Study_Participant_ID) # A slow computation
message(corfit3$consensus.correlation)
g3 <- boxplot(tanh(corfit3$atanh.correlations))

corfit2 <- duplicateCorrelation(logCPM, design=design2, ndups=1, block=meta$Study_Participant_ID) # A slow computation
message(corfit2$consensus.correlation)
g2 <- boxplot(tanh(corfit2$atanh.correlations))

corfit1 <- duplicateCorrelation(logCPM, design=design1, ndups=1, block=meta$Study_Participant_ID) # A slow computation
message(corfit1$consensus.correlation)
g1 <- boxplot(tanh(corfit1$atanh.correlations))

corfit0 <- duplicateCorrelation(logCPM, design=design0, ndups=1, block=meta$Study_Participant_ID) # A slow computation
message(corfit0$consensus.correlation)
g0 <- boxplot(tanh(corfit0$atanh.correlations))

ggsave("dupCor_trainSPA_withBatch.png")

v <- limma::voom(dge, design, plot=TRUE, save.plot = TRUE)
save(v, file="voomOut_trainSPA_withBatch.Rdata")

fit <- lmFit(v, design, correlation = corfit$consensus)
fit <- eBayes(fit)
res = topTable(fit, n=Inf, sort="p", coef=2)

#mapping from ENSEMBL Gene ID to SYMBOL and ENTREZ
ens.str <- substr(rownames(res), 1, 15)
edb <- EnsDb.Hsapiens.v86
res$symbol <- mapIds(edb, keys=ens.str, column="SYMBOL", keytype="GENEID", multiVals="first")
res$entrez <- mapIds(edb, keys=ens.str, column="ENTREZID", keytype="GENEID", multiVals="first")

res.05 <- subset(res, adj.P.Val<0.5) %>% 
  as.data.frame()
write.csv(res.05, file = "results_trainSPA_dupCor_withBatch.05.csv")

#exporting all genes
write.csv(as.data.frame(res), file="results_trainSPA_dupCor_withBatch.all.csv")

#plotting volcano plot
res.05$diffexpressed <- "NO"
res.05$diffexpressed[res.05$logFC > 2 & res.05$adj.P.Val < 0.5] <- "UP"
res.05$diffexpressed[res.05$logFC < -2 & res.05$adj.P.Val < 0.5] <- "DOWN"

g <-  ggplot(data = res, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point(size=0.5) + theme_classic()
#g <- ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point(size=0.5) + theme_classic()
g <- g + theme(panel.grid.major = element_line(color="gray70", size=0.3, linetype=3)) + geom_vline(xintercept=c(-2,2), col="green4") + geom_hline(yintercept=-log10(0.5), col="green4")

#print(g)
ggsave("diff_expressed_trainSPA_dupCor_withBatch.png")

g1 <- g + scale_color_manual(values=c("blue", "gray", "red"))

res.05$delabel <- NA
res.05$delabel[res.05$diffexpressed != "NO"] <- res.05$symbol[res.05$diffexpressed != "NO"]

g1 <- ggplot(data = res.05, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
ggsave("diff_expressed_trainSPA_dupCor_withBatch_label.png")
#print(g1)
