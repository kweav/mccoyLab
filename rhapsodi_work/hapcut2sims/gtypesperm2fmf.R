# Converts GenotypeSperm output for one chromosome to .fmf format for input to hapCUT. Excludes any non-ref or alt observations or observations of multiple alleles from same cell (based on nUMIs, information from GATK table with relevant info for all SNPs)
# By Avery Davis Bell. Use with attribution. (cite preprint/paper)
require(data.table,quietly=T)

##### Functions #####
fread.tfile<-function(tfile){
  # Reads in  VariantsToTable output with chrom, pos, ID, ref, alt; adds column "ind" for index of that SNP
  out<-fread(tfile)
  out[,ind:=1:nrow(out)]
  return(out)
}

tofmf<-function(snptab,daconechr, myn=10){
  # Converts GenotypeSperm output to .fmf-style data.table. Keeps track of observations of SNPs that aren't included (those that observe multiple alleles in one cell, observations of non-ref or alt allele)
  # In: snptab, output of fread.tfile with columns #CHROM, POS, ID, REF, ALT, ind
  #     daconechr, GenotypeSperm output for the one chromosome in snptab; only the following columns: "chr","pos","cell","A","C","G","T","N"
  # Out: list of: fmf, 2 column data.table. First column is cell, second is list of vectors. Each vector is fragment for a consecutive pair of SNPs from that cell. Vector format described in onefrag() comments.
  #               bads, data.table of relevant info of daconechr rows (SNP observations) that were excluded from FMF due to having multiple or the wrong alleles observed
  
  # Subfunctions
  idfails<-function(oneobs){
    # Identifies whether the input SNP observation fails the checks of 1) observing only one allele, 2) seeing only N bases, or 2) seeing a non-ref/non-alt allele
    # In: oneobs, one row representing SNP observation to check. Columns (in order): #CHROM, POS, ID, REF, ALT, ind, chr, cell, A, C, G, T, N
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
  
  onefrag<-function(onepair){
    # Converts one pair of sequential SNP observations to fmf format (exclusive of read ID, which is added later)
    # In: onepair, 2-row datatable of sequential SNP observations. Must have columns POS, REF, ALT, ind, A, C, G, T
    # Out: Either, for SNPs that are not adjacent: vector with - block (2 if SNPs aren't adjacent), off1 (index where first block starts), alleles for SNP in this block in 0/1 alt count format, index of where second block starts, alleles for SNP in this block in 0/1 alt count format, "~~" for fake fastq quality values 
    #      Or, for the less common SNPs that are adjacent:vector with - block (1),off1 (index where first block starts), alleles for SNP in this block in 0/1 alt count format (both), "~~" for fake fastq quality values
    
    
    # format alleles
    if(onepair[1,get(onepair[1,REF])]>0){
      a1<-0
    }else{a1<-1}
    if(onepair[2,get(onepair[2,REF])]>0){
      a2<-0
    }else{a2<-1}
    
    # deal with most cases: SNPs aren't adjacent
    if(onepair[1,ind]<(onepair[2,ind]-1)){
      out<-c(2,onepair[1,ind],a1,onepair[2,ind],a2,"~~")
    }
    
    # deal with adjacent SNPs
    if(onepair[1,ind]==(onepair[2,ind]-1)){
      out<-c(1,onepair[1,ind],paste(a1,a2,sep=""),"~~")
    }
    
    return(out)
  }
  
  onecell.frags<-function(onecell,myn=10){
    # Gets fmf format fragment information for all sequential observations of SNPs in onecell
    # In: onecell, datatable containing observations of all SNPs for one cell. Must have columns POS, REF, ALT, ind, A, C, G, T
    #     myn, number of SNP observations at or below which this cell won't be included. Must be at least 1. For dealing with absence of chromosomes, largely.
    # Out: list with one character vector for each pair. Character vectors structured as described in onefrag() [unless fewer than myn, then null]
    if(nrow(onecell)>myn){  # deal with any cells with only few observations (very occasional aneuploidy)
      out<-lapply(1:(nrow(onecell)-1),function(x) onefrag(onecell[x:(x+1),]))
      return(out)
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
  ckbads<-daconechr[,idfails(.SD),by=list(row.names(daconechr))]$V1 # need to retain as data.table. Unfortunately, this is quite slow. Takes 43 seconds for 10k observations
  bads<-daconechr[ckbads,.(`#CHROM`,POS,ID,REF,ALT,cell,A,C,G,T,N)]
  daconechr<-daconechr[!ckbads,]
  if (sum(ckbads) > 0){ 
    cat(paste("SNP observations to exclude identified (chromosome",daconechr[1,`#CHROM`],")\n")) # for tracking
  }
  
  # Create FMFs on per-cell basis
  setkeyv(daconechr,c("cell","POS"))
  fmf<-daconechr[,list(list(onecell.frags(.SD,myn=myn))),by=cell] # using a list column: vector will display with commas
  fmf<-fmf[sapply(fmf$V1,length)>0,] # get rid of nulls. Now each cell is a list of lists.
  if (nrow(fmf) == 0){
    stop("No fmf conversion possible due to myn argument")
  }
  cat(paste("Initial FMF info generated (chromosome",daconechr[1,`#CHROM`],")\n")) # for tracking
  
  # Return fmf, excludeds
  return(list(fmf,bads))
}

writefmf<-function(fmf,outfile){
  # Disambiguates fmf from 2 column data.table with lists of lists - writes each line including ID based on cell and number
  # In: fmf, second output from tofmf
  #     outfile, path to write this
  
  # Subfunctions
  reformatone<-function(oneout){
    # Converts one row of fmf to a list. Each element in list is correctly-formatted line of output with columns block, fragment ID (cell_#), start position and alleles for appropriate number of blocks, ~~ as quality score
    forout<-unlist(oneout[1,V1],recursive=F) # gets flat list of FMF fragments
    out<-lapply(1:length(forout),function(x) c(forout[[x]][1],paste(oneout$cell,x,sep="_"),forout[[x]][2:length(forout[[x]])])) # adds name to each
    return(out)
  }
  
  # Get flat list with fragment names to output
  towrite<-unlist(sapply(1:nrow(fmf),function(x) reformatone(fmf[x,])),recursive=F)
  
  # write out
  #writeLines(sapply(towrite, paste, collapse="\t"),outfile,sep="\n")
  writeLines(sapply(towrite, paste, collapse=" "), outfile, sep="\n")
}

#### Arguments and input ####
args<-commandArgs(TRUE)
if(length(args)!=4){
  stop("Usage: Rscript gtypesperm2fmf.R
       <1. Path to GATK VariantsToTable output with columns CHROM, POS, ID, REF, ALT. For one chromosome only. >
       <2. Path to GenotypeSperm output - if this file ends in .gz, it will be read as a gzip (via linux; may not work on other platforms depending on zcat install) >
       <3. Path to no-header list of cell barcodes to include in output >
       <4. Output filestem (.fmf file written, as well as GenotypeSperm-format information for SNP observations that were excluded) >
       ",call=F)
}

tfile <- args[1]
dacfile <- args[2]
bcfile <- args[3]
outstem <- args[4]

snptab.glob<-fread.tfile(tfile)
if(substr(dacfile,nchar(dacfile)-1,nchar(dacfile))=="gz"){ # set up command to read gz
  dacfile<-paste('zcat',dacfile,sep=' ')
}
daconechr.glob<-fread(dacfile,select=c("chr","pos","cell","A","C","G","T","N"))[chr==snptab.glob[1,`#CHROM`]&cell%in%fread(bcfile,header=F)$V1] # keep appropriate chromosome, columns of dac only. (this takes time, memory to read in entire DAC and then subset)

#### Convert to .fmf ####
bothouts<-tofmf(snptab.glob,daconechr.glob, myn=1)

# save
write.table(bothouts[[2]],paste(outstem,".excludedsnpobs.txt",sep=""),quote=F,row.names=F,sep="\t")
writefmf(bothouts[[1]],paste(outstem,".fmf",sep=""))