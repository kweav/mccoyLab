library("GenomicAlignments")
library("BiocParallel")
library("GenomicFeatures")
library("Rsamtools")

register(MulticoreParam(20))


gtffile <- file.path("~/resources/WBcel235", "Caenorhabditis_elegans.WBcel235.97.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf")
ebg <- exonsBy(txdb, by = "gene")

## All alignments together

filenamesAll <- list.files(path = "~/workzfs-rmccoy22/kweave23/ykim_rnaseq_20220329/star_out/bam", pattern="*.bam")
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

save(seAll, file = "/scratch/groups/rmccoy22/kweave23/ykim_rnaseq_20220329_summarized_experiment.Rdata")
