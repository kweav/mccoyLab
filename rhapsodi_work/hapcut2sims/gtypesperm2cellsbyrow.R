# Converts GenotypeSperm output for one chromosome to .cellsbyrow.txt format for input to switch computation, CO calling. Excludes any non-ref or alt observations or observations of multiple alleles from same cell (based on nUMIs, information from GATK table with relevant info for all SNPs)
# basically a narrowing. Output format: cell, pos, gt (0 or 1 - count of REF allele)
# By Avery Davis Bell. Use with attribution. (cite preprint/paper)
require(data.table,quietly=T)

##### Functions #####

tocellsbyrow<-function(snptab,daconechr){
  # Converts GenotypeSperm output to .fmf-style data.table. Keeps track of observations of SNPs that aren't included (those that observe multiple alleles in one cell, observations of non-ref or alt allele)
  # In: snptab, output of fread.tfile with columns CHROM, POS, ID, REF, ALT, ind
  #     daconechr, GenotypeSperm output for the one chromosome in snptab; only the following columns: "chr","pos","cell","A","C","G","T","N"
  # Out: list of: fmf, 2 column data.table. First column is cell, second is list of vectors. Each vector is fragment for a consecutive pair of SNPs from that cell. Vector format described in onefrag() comments.
  #               bads, data.table of relevant info of daconechr rows (SNP observations) that were excluded from FMF due to having multiple or the wrong alleles observed
  
  # Subfunctions
  idfails<-function(oneobs){
    # Identifies whether the input SNP observation fails the checks of 1) observing only one allele, 2) seeing only N bases, or 2) seeing a non-ref/non-alt allele
    # In: oneobs, one row representing SNP observation to check. Columns (in order): CHROM, POS, ID, REF, ALT, ind, chr, cell, A, C, G, T, N
    # Out: logical. TRUE if this observation fails either check; FALSE otherwise
    
    if(sum(unlist(oneobs[,.(A,C,G,T)])>0)>1){ # check multiple bases observed
      return(TRUE)
    }
    
    if(oneobs[,N]>0){ # check # N UMIs
      return(TRUE)
    }
    
    if(!(oneobs[,get(oneobs[,REF])]>0|oneobs[,get(oneobs[,ALT])>0])){ # check saw an expected allele
      return(TRUE)
    }
    
    # return if nothing wrong
    return(FALSE)
  }
  
  getgt<-function(oneobs){
    # Returns 0 or 1 - the count of reference allele
    # In: oneobs, one row representing SNP observation to check. Columns (in order): CHROM, POS, ID, REF, ALT, ind, chr, cell, A, C, G, T, N
    # Out: 0 or 1
    if(oneobs[,get(oneobs[,REF])]>0){
      return(0)
    }
    if(oneobs[,get(oneobs[,ALT])]>0){
      return(1)
    }
  }
  
  # Set up so have all relevant data: annotate each observation with SNP table info
  setkey(daconechr,pos)
  setkey(snptab,POS)
  daconechr<-snptab[daconechr]
  if(sum(is.na(daconechr$REF))>0){ # I think I've seen this occasionally
    daconechr<-daconechr[-which(is.na(REF)),]
  }
  
  # Drop any SNP observations of multiple alleles, non-ref or alt alleles
  ckbads<-daconechr[,idfails(.SD),by=row.names(daconechr)]$V1 # need to retain as data.table. Unfortunately, this is quite slow. Takes 43 seconds for 10k observations
  bads<-daconechr[ckbads,.(CHROM,POS,ID,REF,ALT,cell,A,C,G,T,N)]
  daconechr<-daconechr[!ckbads,]
  cat(paste("SNP observations to exclude identified (chromoosome",daconechr[1,CHROM],")\n")) # for tracking
  
  # Get gt
  daconechr[,gt:=sapply(1:nrow(daconechr),function(x) getgt(daconechr[x,]))]
  
  # Get and return clean output
  return(list(daconechr[,.(cell,POS,gt)],bads))
}

#### Arguments and input ####
args<-commandArgs(TRUE)
if(length(args)!=4){
  stop("Usage: Rscript gtypesperm2cellsbyrows.R
       <1. Path to GATK VariantsToTable output with columns CHROM, POS, ID, REF, ALT. For one chromosome only. >
       <2. GenotypeSperm output - if this file ends in .gz, it will be read as a gzip (via linux; may not work on other platforms depending on zcat install) >
       <3. Path to no-header list of cell barcodes to include in output >
       <4. Output filestem (.cellsbyrows.txt file written, as well as GenotypeSperm-format information for SNP observations that were excluded) >
       ",call=F)
}

tfile<-args[1]
dacfile<-args[2]
bcfile<-args[3]
outstem<-args[4]

snptab.glob<-fread(tfile)
if(substr(dacfile,nchar(dacfile)-1,nchar(dacfile))=="gz"){ # set up command to read gz
  dacfile<-paste('zcat',dacfile,sep=' ')
}
daconechr.glob<-fread(dacfile,select=c("chr","pos","cell","A","C","G","T","N"))[chr==snptab.glob[1,CHROM]&cell%in%fread(bcfile,header=F)$V1] # keep appropriate chromosome, columns of dac only. (this takes time, memory to read in entire DAC and then subset)

#### Convert to .fmf ####
bothouts<-tocellsbyrow(snptab.glob,daconechr.glob)

# save
write.table(bothouts[[2]],paste(outstem,".excludedsnpobs.txt",sep=""),quote=F,row.names=F,sep="\t")
setnames(bothouts[[1]],c("cell","pos","gt")) # correct "POS" to be "pos"
setkey(bothouts[[1]],cell) # Sorted by cell
write.table(bothouts[[1]],paste(outstem,".cellsbyrow.txt",sep=""),quote=F,row.names=F,sep="\t")