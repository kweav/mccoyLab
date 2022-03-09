library(SummarizedExperiment)
library(sva)
library(rtracklayer)
library(edgeR)
library(tidyverse)
library(data.table)
library(readxl)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(DESeq2)

load("/home/kweave23/create/batch_correct/create_summarized_experiment_20220207.Rdata")
counts <- as.matrix(assays(seAll)$counts)
for (col in 1:ncol(counts)){
  colnames(counts)[col]<-strsplit(sub("Aligned.sortedByCoord.out.bam", "", colnames(counts)[col]), "_")[[1]][1]
}


#subset to genes with matching gencode ensembl gene IDs on chr1-22
file_gencode <- "/home/kweave23/create/resources/gencode.v34.annotation.gtf"
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

#Importing metadata to differentiate between pregnant and not pregnant and other main covariates
meta <- read.csv("/home/kweave23/create/batch_correct/meta_withSet_kw_20211217.csv", row.names = 1) %>% 
  as.data.frame() %>% 
  mutate(across(c("cDNA_RT_Date", "Library_Prep_Date", "Sequencing_Date", "Study_Participant_ID", "Pregnant"), as.factor))
#subsetting meta data to just the samples with sequencing data            
remove <- setdiff(rownames(meta), colnames(counts))
if (length(which(rownames(meta) %in%  remove)) > 0){
  meta <- meta[-which(rownames(meta) %in%  remove),]
}

#subsetting meta and counts to just training samples
meta <- meta[meta$set == 'train',]
counts <- counts[,which(colnames(counts) %in% rownames(meta))]

#need to sort so name order is the same for setting up design for DESeq experiment
counts_order <- order(colnames(counts))
meta_order <- order(rownames(meta))
counts <- counts[, counts_order]
meta <- meta[meta_order, ]

load("/home/kweave23/create/batch_correct/combatseq_corrected_withshrinkage.RData")


corrected_jg_rts <- corrected_jg_rts[,which(colnames(corrected_jg_rts) %in% rownames(meta))]

dds <- DESeqDataSetFromMatrix(countData = corrected_jg_rts,
                              colData = meta,
                              design = ~ factor(Pregnant) +  factor(cDNA_RT_Date))


keep <- rowSums(counts(dds) >= 10) >= 2
dds<-dds[keep,]

dds  <- DESeq(dds, parallel= TRUE)
res <- results(dds)
filt_threshold <- metadata(res)$filterThreshold
message(filt_threshold)
save(res, file="/home/kweave23/create/batch_correct/res_obj.Rdata")