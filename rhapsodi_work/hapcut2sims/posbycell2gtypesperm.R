# Converts Generative Simulation output to GenotypeSperm output for input to gtypesperm2fmf.R
# GenotypeSperm output will be chr, pos, cell, A, C, G, T, N
# Makes GATK VariantsToTable output with columns CHROM, POS, ID, REF, ALT for input to gtypesperm2fmf.R
# Makes no-header list of cell barcodes to include in output
require(data.table, quietly=T)
require(tidyverse, quietly=T)

##### Arguments and input ####
args <- commandArgs(trailingOnly=TRUE)
pbcfile <- args[1]
outstem <- args[2]


##### Functions #####
fread.mfile <- function(mfile, hetOnly = TRUE){
  # Reads in matrix with first column positions, second column on are different gametes
  # Returns the positions separated from the gamete matrix
  out <- fread(mfile)
  positions <- out[, 1]
  out <- out[,-1]
  if (hetOnly){
    hetRows <- ((rowSums(out == 0, na.rm = TRUE) > 0) & (rowSums(out == 1, na.rm=TRUE) > 0))
    if (sum(hetRows) == 0){
      stop("no hetSNPs so conversion is exiting")
    }
    positions <- positions[hetRows]
    out <- out[hetRows,]
  }
  return(list(positions = positions, out_mat = out))
}

togs <- function(genotypes, barcodes, position, chr = "chrS"){
  nucs <- c("A", "C", "G", "T")
  cells <- c()
  genotypes_list <- list("A" = c(), "C" = c(), "G" = c(), "T" = c(), "N" = c())
  ref <- sample(nucs, size=1)
  alt <- sample(setdiff(nucs, ref), size=1)
  gatkvar_toreturn <- data.frame(`#CHROM`=chr, "POS"=position, "ID" = paste0(chr,"_",position), "REF"=ref, "ALT"=alt, "QUAL" = ".", "FILTER" = "PASS", "INFO"=".", "FORMAT"="GT:GQ", "SAMPLENAME" = "0/1:60", check.names = FALSE)
  for (i in 1:length(barcodes)){
    geno <- genotypes[i]
    if (!is.na(geno)){
      cells <- c(cells, barcodes[i])
      if (geno == 1){
        genotypes_list[alt][[1]] <- c(genotypes_list[alt][[1]],1)
        for (nuc in setdiff(nucs,alt)){
          genotypes_list[nuc][[1]] <- c(genotypes_list[nuc][[1]],0)
        }
      } else if (geno == 0){
        genotypes_list[ref][[1]] <- c(genotypes_list[ref][[1]],1)
        for (nuc in setdiff(nucs,ref)){
          genotypes_list[nuc][[1]] <- c(genotypes_list[nuc][[1]],0)
        }
      }
      genotypes_list["N"][[1]] <- c(genotypes_list["N"][[1]],0)
    }
  }
  if (length(cells) > 0){
    to_return <- data.frame("chr"=chr, "pos"=position, "cell"=cells, "A" = genotypes_list["A"][[1]], "C" = genotypes_list["C"][[1]], "G" = genotypes_list["G"][[1]], "T" = genotypes_list["T"[[1]]], "N" = genotypes_list["N"][[1]])
  } else {
    to_return <- data.frame()
  }
  return(list(gs = to_return, gatkvar = gatkvar_toreturn))
  #why does gs only have 7 rows?
}

mattogs <- function(mat, positions){
  gs <- data.frame()
  gatkvar <- data.frame()
  for (x in 1:length(positions)){
    out <- togs(unlist(mat[x,]), colnames(mat), positions[x][[1]])
    gs <- rbind(gs, out$gs )
    gatkvar <- rbind(gatkvar, out$gatkvar)
  }
  #rbindlist two separate dataframes in a returned list?s
  return(list(gs = gs, gatkvar = gatkvar))
}

set.seed(789)
#read_out <- fread.mfile("runGen_gam_1000_snp_30000_cov_0.001_seqerr_0.005_avgr_1_rs_42_gametedf_na_truth_afseqednm.csv")
read_out <- fread.mfile(pbcfile)
togs_out <- mattogs(read_out$out_mat, read_out$positions[[1]])
cell_barcodes <- colnames(read_out$out_mat)

#write.table(togs_out$gs, file = "~/mccoyLab/rhapsodi_work/hapcut2sims/gtypesperm_g1000_s30000_cov_0.001_seqerr_0.005_avgr_1_rs42.txt", row.names = FALSE)
write.table(togs_out$gs, file = paste0(outstem, "_gtypesperm.txt"), row.names = FALSE)
#write.table(togs_out$gatkvar, file="~/mccoyLab/rhapsodi_work/hapcut2sims/variantinfo_g1000_s30000_cov_0.001_seqerr_0.005_avgr_1_rs42.vcf", row.names = FALSE, sep = "\t", quote=FALSE)
write.table(togs_out$gatkvar, file = paste0(outstem, "_variantinfo.vcf"), row.names = FALSE, sep = "\t", quote=FALSE)
#write.table(cell_barcodes, file="~/mccoyLab/rhapsodi_work/hapcut2sims/cellbarcodes_g1000_s30000_cov_0.001_seqerr_0.005_avgr_1_rs42.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(cell_barcodes, file=paste0(outstem, "_cellbarcodes.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)