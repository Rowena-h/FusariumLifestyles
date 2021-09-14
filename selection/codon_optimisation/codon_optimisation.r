#!/usr/bin/env Rscript
##Script to parse OrthoFinder and ribosomal protein BLAST results to calculate codon optimisation for core single copy genes##

library(seqinr)
library(coRdon)
library(tAI)

args <- commandArgs(trailingOnly=TRUE)

#Test if there is one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied: the OrthoFinder results directory (ending in a forward slash)", call.=FALSE)
} 

dir <- args[1]


##GC12 GC3 content for neutrality plot##

#List all core single-copy codon alignments
codon.files <- list.files("../alignments/codon/", pattern=".fa")

#Make empty list for GC content results
gc.list <- list()

#For each core SC gene...
for (i in codon.files) {
  
  #Try to read in codon alignment
  alignment <- tryCatch(read.fasta(paste0("../alignments/codon/", i)), error=function(e) NULL)
  
  #If the alignment is not empty...
  if (lengths(alignment[1]) > 0) {
    #For each taxon in the alignment...
    for (j in names(alignment)) {
      
      #Calculate GC12
      gc.list[[j]][[i]]["GC12"] <- (GC1(alignment[[j]])  + GC2(alignment[[j]])) / 2
      #Calculate GC3
      gc.list[[j]][[i]]["GC3"] <- GC3(alignment[[j]])
      
    }
  }
  
}

#Convert to dataframe
gc.df <- data.frame(t(as.data.frame(gc.list)))
gc.df$taxon <- sub(".OG000.*", "", rownames(gc.df))


##Codon optimisation##

#Read in orthogroups from OrthoFinder
orthogroups <- read.csv(paste0(dir, "Orthogroups/Orthogroups.tsv"), row.names=1, sep="\t", check.names=FALSE)

#Read in 'unassigned genes' i.e. species specific genes and combines with orthogroups dataframe
unassigned <- read.csv(paste0(dir, "Orthogroups/Orthogroups_UnassignedGenes.tsv"),
                       row.names=1, sep="\t", check.names=FALSE)

orthogroups <- rbind(orthogroups, unassigned)

#For each taxon...
print("Reading in ribosomal proteins")
for (i in colnames(orthogroups)) {
  
  #Read in the list of ribosomal proteins
  ribosomes <- scan(paste0(i, ".faa_blast_ribosomes"), character(), quote="")
  #Replace pipes (|) with hyphens
  ribosomes <- gsub("\\|", "-", ribosomes)
  assign(paste0(i, ".ribosomes"), ribosomes)
  
  #Replace pipes (|) with hyphens
  orthogroups[,i] <- gsub("\\|", "-", orthogroups[,i])
  
}


#Make dataframe with gene counts for each orthogroup
orthogroups.copies <- orthogroups

for (i in 1:length(colnames(orthogroups.copies))) {
  orthogroups.copies[, i] <- sapply(strsplit(orthogroups[, i], " "), length)
}

#Make dataframe for ribosomal protein counts
ribosome.count <- data.frame(matrix(0, ncol=ncol(orthogroups), nrow=nrow(orthogroups)))
colnames(ribosome.count) <- colnames(orthogroups)
rownames(ribosome.count) <- rownames(orthogroups)

#For each taxon...
for (i in 1:length(colnames(ribosome.count))) {
  
  #Print progress
  cat("Counting number of ribosomal proteins in each orthogroup: ", (i - 1), "/", length(colnames(ribosome.count)), " taxa", "\r")
  
  #Retrieve the list of ribosomal proteins
  ribosomes <- get(paste0(colnames(ribosome.count)[i], ".ribosomes"))
  
  #For each row in the list (i.e. protein)...
  for (j in 1:length(ribosomes)) {
    
    #Retrieve orthogroup ID
    ribosome <- grep(ribosomes[j], orthogroups[, i])
    
    #Search for orthogroup ID in corresponding column of orthogroups dataframe and add 1 to orthogroup count
    ribosome.count[ribosome, i] <- ribosome.count[ribosome, i] + 1
    
  }
  
}
print(paste0("Counting number of ribosomal proteins in each orthogroup: ", i, "/", length(colnames(ribosome.count)), " taxa"))


#For each taxon...
for (i in colnames(orthogroups.copies)) {
  
  #Print progress
  cat("Generating codon tables and RSCU values ", (which(colnames(orthogroups.copies)  == i) - 1),
      "/", length(colnames(orthogroups.copies)), " taxa", "\r")
  
  #Read in core single-copy proteins
  prots <- readSet(file=paste0(i, "_coreSC.fa"))
  #Make table of codon counts for core SC proteins
  codon.table <- codonTable(prots)
  #Calculate relative synonymous codon usage
  rscu <- uco(unlist(strsplit(paste(as.vector(prots), collapse=""), "")), index="rscu")
    
  assign(paste0("rscu.", i), rscu)
  assign(paste0("codon.table.", i), codon.table)
  
}
print(paste0("Generating codon tables and RSCU values ", which(colnames(orthogroups.copies)  == i), "/", length(colnames(ribosome.count)), " taxa"))

#Make empty vector to label which core SC proteins are ribosomal
test.set <- rep(FALSE, length(codon.table))
ribosome.orthos <- rownames(ribosome.count)[which(rowSums(ribosome.count) > 0)]
test.set[na.omit(match(ribosome.orthos, names(prots)[which(lengths(prots) > 0)]))] <- TRUE

#Read in orthogroup data
load("../../CSEP_prediction/orthogroup-matrices-2021-09-09.RData")

#Core, single-copy CSEPs
core.SC.mixed <- Reduce(intersect,
                        list(orthogroups.stats.ingroup0$orthogroup[which(
                          orthogroups.stats.ingroup0$copy_number == "single")],
                          orthogroups.stats.ingroup0$orthogroup[which(
                            orthogroups.stats.ingroup0$category == "core")],
                          orthogroups.stats.ingroup0$orthogroup[which(
                            orthogroups.stats.ingroup0$CSEP != "")]))

#Make empty vector for GC3 content
gc3.list <- list()

#Make empty dataframe for codon optimisation (S) results
s.df <- data.frame(taxon=colnames(orthogroups.copies),
                   S=NA,
                   S.CSEP=NA,
                   S.other=NA)

cai.list <- list()
enc.list <- list()

#For each taxon...
for (i in colnames(orthogroups.copies)) {
  
  #Print progress
  cat("Calculating codon statistics for each core SC orthogroup: ",
      (which(colnames(orthogroups.copies) == i) - 1), "/",
      length(colnames(orthogroups.copies)), " taxa", "\r")
  
  codon.table <- get(paste0("codon.table.", i))
  
  #Calculate codon adaptation index
  cai <- CAI(codon.table, subsets=list(ribosomes=test.set), stop.rm=TRUE)
  
  cai.list[[i]] <- cai
  
  #Calculate effective number of codons
  enc <- ENC(codon.table)
  
  enc.list[[i]] <- enc
  
  #Read in core SC sequences
  fasta <- tryCatch(read.fasta(file=paste0(i, "_coreSC.fa")), error=function(e) NULL)
  
  #For each core SC orthogroup...
  for (j in names(fasta)) {
    if (length(fasta[[j]]) > 0) {
      
      #Calculate GC3 content
      gc3.list[[i]][j] <- GC3(fasta[[j]])
      
    }
    
  }
  
  #Calculate S
  s.df$S[s.df$taxon == i] <- get.s(cai, enc, as.vector(gc3.list[[i]]))
  #...for CSEPs
  s.df$S.CSEP[s.df$taxon == i] <-get.s(cai[match(core.SC.mixed, getID(codon.table))],
                                           enc[match(core.SC.mixed, getID(codon.table))],
                                           as.vector(gc3.list[[i]])[match(core.SC.mixed, getID(codon.table))])
  #...for non-CSEPs
  s.df$S.other[s.df$taxon == i] <-get.s(cai[-match(core.SC.mixed, getID(codon.table))],
                                           enc[-match(core.SC.mixed, getID(codon.table))],
                                           as.vector(gc3.list[[i]])[-match(core.SC.mixed, getID(codon.table))])
  
}
print(paste0("Calculating codon statistics for each core SC orthogroup: ",
             which(colnames(orthogroups.copies) == i), "/",
             length(colnames(orthogroups.copies)), " taxa"))

print("Testing for significant difference in S values between CSEPs and other genes:")

shapiro.e <- shapiro.test(s.df$S.CSEP)
shapiro.o <- shapiro.test(s.df$S.other)

if (shapiro.e$p.value < 0.05 || shapiro.o$p.value < 0.05) {
  print("Reject normality, doing Wilcoxon test")
  wilcox.test(x=s.df$S.CSEP, y=s.df$S.other)
} else {
  print("Normal data, doing t-test")
  t.test(x=s.df$S.CSEP, y=s.df$S.other)
}

print(paste0("Results saved in codon_optimisation-", Sys.Date(), ".RData"))
save(list=c(ls(pattern="rscu\\."), "s.df", "gc.df"), file=paste0("codon_optimisation-", Sys.Date(), ".RData")) 
