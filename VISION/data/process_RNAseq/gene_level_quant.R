library(here)
library(tximport)
library(rhdf5)
library(rlist)

#Usage: (r_env) kweave23@comp2:/project/vision/Data/RNAseq/rnaseq_quant$ date;time Rscript ../gene_level_quant.R &> ../quant_gene_level.txt

dir <- here()
tx2gene <- read.csv('../transcript_to_gene_df.Mus_musculus.GRCm38.96.csv', header=FALSE)
files_list <- list.files(dir)
files <- file.path(dir, files_list, "abundance.h5")
names(files) <- files_list
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene, txOut = FALSE) #this will summarize at the gene-level
list.save(txi.kallisto, '../txi.kallisto.list.rdata')
#txi.kallisto <- list.load('../txi.kallisto.list.rdata')
#typeof(txi.kallisto$abundance)
write.csv(txi.kallisto$abundance, '../txi.kallisto.abundances.csv')