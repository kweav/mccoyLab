library("GenomicAlignments")
library("BiocParallel")
library("GenomicFeatures")
library("Rsamtools")


register(MulticoreParam(20))


gtffile <- file.path("/scratch/groups/rmccoy22/kweave23/resources", "gencode.v34.annotation.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf")
ebg <- exonsBy(txdb, by = "gene")

## All alignments together

filenamesAll <- list.files(path = "/scratch/groups/rmccoy22/kweave23/create_outcome_pred/alignments", pattern="*.bam")
fileAllexists <- unlist(lapply(1:length(filenamesAll), function(x) file.exists(filenamesAll[x])))
message(sum(fileAllexists))
message(length(filenamesAll))
bamfilesAll <- BamFileList(filenamesAll, yieldSize = 2000000)
length(bamfilesAll)

seAll <- summarizeOverlaps(features = ebg, 
                           reads = bamfilesAll,
                           mode = "Union",
                           singleEnd = FALSE,
                           ignore.strand = TRUE,
                           fragments = TRUE)

save(seAll, file = "/scratch/groups/rmccoy22/kweave23/create_outcome_pred/create_summarized_experiment_allNov2021.Rdata")