#!/usr/bin/env Rscript
##Script to parse OrthoFinder and CSEPfilter results into presence-absence matrix of orthogroups for each taxon##

library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

#Test if there is one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied: the OrthoFinder results directory (ending in a forward slash)", call.=FALSE)
} 

dir <- args[1]

#Read in orthogroups from OrthoFinder
orthogroups <- read.csv(paste0(dir, "Orthogroups/Orthogroups.tsv"), row.names=1, sep="\t", check.names=FALSE)

#Read in 'unassigned genes' - i.e. strain specific genes - and combine with orthogroups dataframe
unassigned <- read.csv(paste0(dir, "Orthogroups/Orthogroups_UnassignedGenes.tsv"),
                       row.names=1, sep="\t", check.names=FALSE)

orthogroups <- rbind(orthogroups, unassigned)

#For each sample...
message("Reading in CSEPs")
for (i in colnames(orthogroups)) {
  
  #Read in the list of CSEPs
  CSEPs <- scan(paste0(i, ".faa_candidate_effectors"), character(), quote="")
  #Replace pipes (|) with hyphens
  CSEPs <- gsub("\\|", "-", CSEPs)
  orthogroups[,i] <- gsub("\\|", "-", orthogroups[,i])
  
  assign(paste0(i, ".CSEPs"), CSEPs)
  
}

#Make dataframe with gene counts for each orthogroup
orthogroups.copies <- orthogroups

for (i in 1:length(colnames(orthogroups.copies))) {
  orthogroups.copies[, i] <- sapply(strsplit(orthogroups[, i], " "), length)
}

#Make dataframe for CSEP counts
CSEP.count <- data.frame(matrix(0, ncol=ncol(orthogroups), nrow=nrow(orthogroups)))
colnames(CSEP.count) <- colnames(orthogroups)
rownames(CSEP.count) <- rownames(orthogroups)

message("Counting number of CSEPs in each orthogroup:")

#For each sample...
for (i in 1:length(colnames(CSEP.count))) {
  
  #Print progress
  message((i - 1), "/", length(colnames(CSEP.count)))
  
  #Retrieve the list of potential CSEPs
  CSEPs <- get(paste0(colnames(CSEP.count)[i], ".CSEPs"))
  
  #For each CSEP...
  for (j in 1:length(CSEPs)) {
    
    #Retrieve CSEP
    CSEP <- grep(CSEPs[j], orthogroups[, i])
    
    #Search for CSEP in corresponding column of orthogroups dataframe and add 1 to orthogroup count
    CSEP.count[CSEP, i] <- CSEP.count[CSEP, i] + 1
    
  }
  
}
message(i, "/", length(colnames(CSEP.count)))

#Read in sample metadata
metadata <- read.csv("../metadata.csv")

#Including and excluding the outgroup...
for (ingroup in c(0, 1)) {
  
  if (ingroup == 1) {
    
    orthogroups.copies.ingroup <- orthogroups.copies[which(!is.na(match(colnames(orthogroups.copies), metadata$file2[metadata$ingroup == ingroup])))]
    CSEP.count.ingroup <- CSEP.count[which(!is.na(match(colnames(CSEP.count), metadata$file2[metadata$ingroup == ingroup])))]
    
  } else {
    
    orthogroups.copies.ingroup <- orthogroups.copies
    CSEP.count.ingroup <- CSEP.count
    
  }
  
  #Make dataframe to collect stats for orthogroups
  orthogroups.stats <- data.frame(orthogroup=rownames(orthogroups.copies.ingroup),
                                  copy_number="multi",
                                  category=NA,
                                  CSEP=NA,
                                  stringsAsFactors=FALSE)
  
  #Add whether orthogroups are single copy
  orthogroups.stats$copy_number[match(rownames(orthogroups.copies.ingroup[apply(orthogroups.copies.ingroup < 2, 1, all),]), orthogroups.stats$orthogroup)] <- "single"
  
  #Categorise orthogroups according to sharedness
  message("Assigning orthogroups as core, accessory or specific for ingroup ", ingroup, ":")
  #Create bar to show progress
  progress.bar <- txtProgressBar(1, length(rownames(orthogroups.copies.ingroup)), initial=0, char="=", style=3)
  for (j in 1:length(rownames(orthogroups.copies.ingroup))) {
    
    #Update progress bar
    setTxtProgressBar(progress.bar, j)
    
    if (length(which(orthogroups.copies.ingroup[j,] == 0)) == (length(colnames(orthogroups.copies.ingroup)) - 1)) {
      orthogroups.stats$category[orthogroups.stats$orthogroup == rownames(orthogroups.copies.ingroup)[j]] <-
        "specific"
    }
    if (length(which(orthogroups.copies.ingroup[j,] == 0)) == 0) {
      orthogroups.stats$category[orthogroups.stats$orthogroup == rownames(orthogroups.copies.ingroup)[j]] <-
        "core"
    }
    if (length(which(orthogroups.copies.ingroup[j,] == 0)) < (length(colnames(orthogroups.copies.ingroup)) - 1) && length(which(orthogroups.copies.ingroup[j,] == 0)) > 0) {
      orthogroups.stats$category[orthogroups.stats$orthogroup == rownames(orthogroups.copies.ingroup)[j]] <-
        "accessory"
    }
  }
  close(progress.bar)
  
  #Determine which orthogroups contain only CSEPs (SP-only) or both CSEPs and other genes (SP-mixed) 
  CSEPs <- names(which(rowSums(CSEP.count.ingroup) > 0))
  
  CSEP.type <- list()
  
  message("Assigning CSEP orthogroups as CSEP-only or CSEP-mixed for ingroup ", ingroup, ":")
  
  #For each sample...
  for (i in 1:length(colnames(orthogroups.copies.ingroup))) {
    
    #Print progress
    message((i - 1), "/", length(colnames(orthogroups.copies.ingroup)))
    
    for (j in CSEPs) {
      if (orthogroups.copies.ingroup[j,i] == CSEP.count.ingroup[j,i]) {
        CSEP.type[[j]][i] <- "CSEP"
      } else {
        CSEP.type[[j]][i] <- "mixed"
      }
    }
    
  }
  message(i, "/", length(colnames(orthogroups.copies.ingroup)))
  
  for (j in CSEPs) {
    if(grepl("mixed", CSEP.type[j]) == TRUE) {
      orthogroups.stats$CSEP[orthogroups.stats$orthogroup == j] <- "mixed"
    } else {
      orthogroups.stats$CSEP[orthogroups.stats$orthogroup == j] <- "CSEP"
    }
  }
  
  assign(paste0("orthogroups.stats.ingroup", ingroup), orthogroups.stats)
  assign(paste0("orthogroups.copies.ingroup", ingroup), orthogroups.copies.ingroup)
  assign(paste0("CSEP.count.ingroup", ingroup), CSEP.count.ingroup)
  
}

message(paste0("Results saved in orthogroup-matrices-", Sys.Date(), ".RData"))
save(orthogroups.stats.ingroup0,
     orthogroups.stats.ingroup1,
     orthogroups.copies.ingroup0,
     orthogroups.copies.ingroup1,
     CSEP.count.ingroup0,
     CSEP.count.ingroup1,
     file=paste0("orthogroup-matrices-", Sys.Date(), ".RData"))
