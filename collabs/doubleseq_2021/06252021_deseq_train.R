library(DESeq2)
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

load("/home/kweave23/create/create_summarized_experiment.Rdata")
counts <- as.matrix(assays(se)$counts)
for (col in 1:ncol(counts)){
  colnames(counts)[col]<-strsplit(sub("Aligned.sortedByCoord.out.bam", "", colnames(counts)[col]), "_")[[1]][1]
}

#subset to just training data
embryoID_by_set <- read_csv("/home/kweave23/create/embryoID_by_set.csv")
train_embryoIDs = embryoID_by_set$embryoID[embryoID_by_set$set == "train"]
train_cols <- which(colnames(counts) %in% train_embryoIDs)

counts <- counts[,train_cols]

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

#Importing metadata to differentiate between pregnant and not pregnant
meta <- read.csv("/home/kweave23/create/tidied_meta_CREATE_kw.csv", row.names = 1) %>% 
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

#Filtering out counts 0 and 1
keep <- rowSums(countdata) > 1
dds <- countdata[keep, ]

dds2 <- DESeqDataSetFromMatrix(countData = dds, colData = coldata, design = ~ Pregnant + Sequencing_Date)
#apply a variance stabilizing transformation
vsd <- vst(dds2, blind = FALSE)

#Filtering out counts 0 and 1
keep <- rowSums(counts(dds2)) > 1
dds2 <- dds2[keep, ]

#Running differential expression
dds2 <- DESeq(dds2)
res <- results(dds2)

#mapping from ENSEMBL Gene ID to SYMBOL and ENTREZ
ens.str <- substr(rownames(res), 1, 15)
edb <- EnsDb.Hsapiens.v86
res$symbol <- mapIds(edb, keys=ens.str, column="SYMBOL", keytype="GENEID", multiVals="first")
res$entrez <- mapIds(edb, keys=ens.str, column="ENTREZID", keytype="GENEID", multiVals="first")

#exporting significant genes
res.05 <- subset(res, padj<0.05) %>% 
  as.data.frame()
write.csv(res.05, file = "results_trainSPA_withBatch.05.csv")

#exporting all genes
write.csv(as.data.frame(res), file="results_trainSPA_withBatch.all.csv")

#plotting volcano plot
res.05$diffexpressed <- "NO"
res.05$diffexpressed[res.05$log2FoldChange > 2 & res.05$padj < 0.05] <- "UP"
res.05$diffexpressed[res.05$log2FoldChange < -2 & res.05$padj < 0.05] <- "DOWN"

g <- ggplot(data=res.05, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + geom_point(size=0.5) + theme_classic()
g <- g + theme(panel.grid.major = element_line(color="gray70", size=0.3, linetype=3)) + geom_vline(xintercept=c(-2,2), col="green4") + geom_hline(yintercept=-log10(0.05), col="green4")
#print(g)
ggsave("diff_expressed_trainSPA_with_batch.png")

g1 <- g + scale_color_manual(values=c("blue", "gray", "red"))

res.05$delabel <- NA
res.05$delabel[res.05$diffexpressed != "NO"] <- res.05$symbol[res.05$diffexpressed != "NO"]

g1 <- ggplot(data = res.05, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
ggsave("diff_expressed_trainSPA_with_batch_label.png")
#print(g1)