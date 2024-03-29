#!/usr/bin/env Rscript
##Script to parse OrthoFinder and ribosomal protein BLAST results to calculate codon optimisation for core single copy genes##

library(seqinr)
library(coRdon)
library(tAI)
library(rstatix)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

#Test if there is one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied: the OrthoFinder results directory (ending in a forward slash)", call.=FALSE)
} 

dir <- args[1]


##Codon optimisation##

#Read in orthogroups from OrthoFinder
orthogroups <- read.csv(paste0(dir, "Orthogroups/Orthogroups.tsv"), row.names=1, sep="\t", check.names=FALSE)

#Read in 'unassigned genes' i.e. species specific genes and combines with orthogroups dataframe
unassigned <- read.csv(paste0(dir, "Orthogroups/Orthogroups_UnassignedGenes.tsv"),
                       row.names=1, sep="\t", check.names=FALSE)

orthogroups <- rbind(orthogroups, unassigned)

#For each taxon...
message("Reading in ribosomal proteins")
for (i in colnames(orthogroups)) {
  
  #Read in the list of ribosomal proteins
  ribosomes <- scan(paste0(i, ".faa_ribosomes"), character(), quote="")
  #Replace pipes (|) with hyphens
  ribosomes <- gsub("\\|", "-", ribosomes)
  assign(paste0(i, ".ribosomes"), ribosomes)
  
  #Replace pipes (|) with hyphens
  orthogroups[,i] <- gsub("\\|", "-", orthogroups[,i])
  
}


#Make dataframe for ribosomal protein counts
ribosome.count <- data.frame(matrix(0, ncol=ncol(orthogroups), nrow=nrow(orthogroups)))
colnames(ribosome.count) <- colnames(orthogroups)
rownames(ribosome.count) <- rownames(orthogroups)

message("Counting number of ribosomal proteins in each orthogroup:")

#For each taxon...
for (i in 1:length(colnames(ribosome.count))) {
  
  #Print progress
  message((i - 1), "/", length(colnames(ribosome.count)))
  
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
message(i, "/", length(colnames(ribosome.count)))

message("Generating codon tables and RSCU values:")

#For each taxon...
for (i in colnames(orthogroups)) {
  
  #Print progress
  message((which(colnames(orthogroups)  == i) - 1),
          "/", length(colnames(orthogroups)))
  
  #Read in core single-copy proteins
  prots <- readSet(file=paste0(i, "_coreSC.fa"))
  #Make table of codon counts for core SC proteins
  codon.table <- codonTable(prots)
  #Calculate relative synonymous codon usage
  rscu <- uco(unlist(strsplit(paste(as.vector(prots), collapse=""), "")), index="rscu")
  
  assign(paste0("rscu.", i), rscu)
  assign(paste0("codon.table.", i), codon.table)
  
}
message((which(colnames(orthogroups)  == i)),
        "/", length(colnames(orthogroups)))

#Make empty vector to label which core SC proteins are ribosomal
test.set <- rep(FALSE, length(codon.table))
ribosome.orthos <- rownames(ribosome.count)[which(rowSums(ribosome.count) > 0)]
test.set[na.omit(match(ribosome.orthos, names(prots)[which(lengths(prots) > 0)]))] <- TRUE

#Read in orthogroup data
load("../../CSEP_CAZyme_prediction/orthogroup-matrices-2022-02-10.RData")

#Core, single-copy CSEPs
core.SC.csepmixed <- Reduce(intersect,
                            list(orthogroups.stats.ingroup0$orthogroup[which(
                              orthogroups.stats.ingroup0$copy_number == "single")],
                              orthogroups.stats.ingroup0$orthogroup[which(
                                orthogroups.stats.ingroup0$category == "core")],
                              orthogroups.stats.ingroup0$orthogroup[which(
                                !is.na(orthogroups.stats.ingroup0$CSEP))]))

#Core, single-copy CAZymes
core.SC.cazymemixed <- Reduce(intersect,
                              list(orthogroups.stats.ingroup0$orthogroup[which(
                                orthogroups.stats.ingroup0$copy_number == "single")],
                                orthogroups.stats.ingroup0$orthogroup[which(
                                  orthogroups.stats.ingroup0$category == "core")],
                                orthogroups.stats.ingroup0$orthogroup[which(
                                  !is.na(orthogroups.stats.ingroup0$CAZyme))]))

#Make empty vector for GC3 content
gc3.list <- list()

#Make empty dataframe for codon optimisation (S) results
s.df <- data.frame(taxon=colnames(orthogroups),
                   S=NA,
                   S.CSEP=NA,
                   S.CAZyme=NA,
                   S.other=NA)

cai.list <- list()
enc.list <- list()

message("Calculating codon statistics for each core SC orthogroup:")

#For each taxon...
for (i in colnames(orthogroups)) {
  
  #Print progress
  message((which(colnames(orthogroups) == i) - 1), "/",
          length(colnames(orthogroups)))
  
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
  s.df$S.CSEP[s.df$taxon == i] <-
    get.s(cai[match(core.SC.csepmixed, getID(codon.table))],
          enc[match(core.SC.csepmixed, getID(codon.table))],
          as.vector(gc3.list[[i]])[match(core.SC.csepmixed,
                                         getID(codon.table))])
  #...for CAZymes
  s.df$S.CAZyme[s.df$taxon == i] <-
    get.s(cai[match(core.SC.cazymemixed, getID(codon.table))],
          enc[match(core.SC.cazymemixed, getID(codon.table))],
          as.vector(gc3.list[[i]])[match(core.SC.cazymemixed,
                                         getID(codon.table))])
  #...for non-CSEPs/CAZymes
  s.df$S.other[s.df$taxon == i] <-
    get.s(cai[-match(union(core.SC.csepmixed, core.SC.cazymemixed),
                     getID(codon.table))],
          enc[-match(union(core.SC.csepmixed, core.SC.cazymemixed),
                     getID(codon.table))],
          as.vector(gc3.list[[i]])[-match(union(core.SC.csepmixed, core.SC.cazymemixed),
                                          getID(codon.table))])
  
}
message((which(colnames(orthogroups) == i)), "/",
        length(colnames(orthogroups)))

message(paste0("Results saved in codon_optimisation-", Sys.Date(), ".RData"))
save(list=c(ls(pattern="rscu\\."), "s.df"), file=paste0("codon_optimisation-", Sys.Date(), ".RData")) 
