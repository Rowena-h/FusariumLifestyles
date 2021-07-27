#!/usr/bin/env Rscript
##Script to parse OrthoFinder and ribosomal protein BLAST results to calculate codon optimisation for core single-copy genes##

library(seqinr)
library(coRdon)
library(tAI)

args <- commandArgs(trailingOnly=TRUE)

#Test if there is one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied: the OrthoFinder results directory (ending in a forward slash)", call.=FALSE)
} 

dir <- args[1]

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
  cat("Counting number of effectors in each orthogroup: ", (i - 1), "/", length(colnames(ribosome.count)), " taxa", "\r")
  
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
  
  #Read in core single-copy proteins
  prots <- readSet(file=paste0(i, "_coreSC.fa"))
  #Make table of codon counts for core SC proteins
  codon.table <- codonTable(prots)
    
  assign(paste0("codon.table.", i), codon.table)
  
}

#Make empty vector to label which core SC proteins are ribosomal
test.set <- rep(FALSE, length(codon.table))
ribosome.orthos <- rownames(ribosome.count)[which(rowSums(ribosome.count) > 0)]
test.set[na.omit(match(ribosome.orthos, names(prots)[which(lengths(prots) > 0)]))] <- TRUE

#Make empty vector for GC3 content
gc3.list <- list()

#Make empty dataframe for codon optimisation (S) results
s.df <- data.frame(taxon=colnames(orthogroups.copies),
                   S=NA)

#For each taxon...
for (i in colnames(orthogroups.copies)) {
  
  #Print progress
  cat("Calculating codon statistics for each core SC orthogroup: ",
      (which(colnames(orthogroups.copies) == i) - 1), "/",
      length(colnames(orthogroups.copies)), " taxa", "\r")
  
  #Calculate codon adaptation index
  cai <- CAI(get(paste0("codon.table.", i)), subsets=list(ribosomes=test.set), stop.rm=TRUE)
  
  #Calculate effective number of codons
  enc <- ENC(get(paste0("codon.table.", i)))
  
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
  
}
print(paste0("Calculating codon statistics for each core SC orthogroup: ",
             which(colnames(orthogroups.copies) == i), "/",
             length(colnames(orthogroups.copies)), " taxa"))


print(paste0("Results saved in codon_optimisation-", Sys.Date(), ".RData"))
save(s.df, file=paste0("codon_optimisation-", Sys.Date(), ".RData")) 
