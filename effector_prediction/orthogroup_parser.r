#!/usr/bin/env Rscript
##Script to parse OrthoFinder and SPfilter results into a matrix of effector counts for each taxon##

library(dplyr)

args=commandArgs(trailingOnly=TRUE)

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

#For each sample...
print("Reading in candidate effectors")
for (i in colnames(orthogroups)) {
  
  #Read in the list of candidate effectors
  effectors <- scan(paste0(i, ".faa_candidate_effectors"), character(), quote="")
  #Replace pipes (|) with hyphens
  effectors <- gsub("\\|", "-", effectors)
  assign(paste0(i, ".effectors"), effectors)
  
  #Replace pipes (|) with hyphens
  orthogroups[,i] <- gsub("\\|", "-", orthogroups[,i])
  
}

#Make dataframe with gene counts for each orthogroup
orthogroups.copies <- orthogroups

for (i in 1:length(colnames(orthogroups.copies))) {
  orthogroups.copies[, i] <- sapply(strsplit(orthogroups[, i], " "), length)
}

#Make dataframe for effector counts
effector.count <- data.frame(matrix(0, ncol=ncol(orthogroups), nrow=nrow(orthogroups)))
colnames(effector.count) <- colnames(orthogroups)
rownames(effector.count) <- rownames(orthogroups)

#For each sample...
for (i in 1:length(colnames(effector.count))) {
  
  #Print progress
  cat("Counting number of effectors in each orthogroup: ", (i - 1), "/", length(colnames(effector.count)), " taxa", "\r")
  
  #Retrieve the list of potential effectors
  effectors <- get(paste0(colnames(effector.count)[i], ".effectors"))
  
  #For each row in the list (i.e. potential effector)...
  for (j in 1:length(effectors)) {
    
    #Retrieve effector
    effector <- grep(effectors[j], orthogroups[, i])
    
    #Search for effector in corresponding column of orthogroups dataframe and add 1 to orthogroup count
    effector.count[effector, i] <- effector.count[effector, i] + 1
    
  }
  
}
print(paste0("Counting number of effectors in each orthogroup: ", i, "/", length(colnames(effector.count)), " taxa"))

#Read in sample metadata
metadata <- read.csv("../metadata.csv")

for (ingroup in c(0, 1, 2)) {
  
  if (ingroup == 1) {
    
    orthogroups.copies.ingroup <- orthogroups.copies[which(!is.na(match(colnames(orthogroups.copies), metadata$file2[metadata$ingroup == ingroup])))]
    effector.count.ingroup <- effector.count[which(!is.na(match(colnames(effector.count), metadata$file2[metadata$ingroup == ingroup])))]
    
  } else if (ingroup == 2) {
    
    orthogroups.copies.ingroup <- orthogroups.copies[which(!is.na(match(colnames(orthogroups.copies), metadata$file2[metadata$ingroup != "outgroup"])))]
    effector.count.ingroup <- effector.count[which(!is.na(match(colnames(effector.count), metadata$file2[metadata$ingroup != "outgroup"])))]
    
  } else {
    
    orthogroups.copies.ingroup <- orthogroups.copies
    effector.count.ingroup <- effector.count
    
  }

  orthogroups.stats <- data.frame(orthogroup=rownames(orthogroups.copies.ingroup),
                                  copy_number="multi",
                                  secretome=NA,
                                  effector=NA,
                                  stringsAsFactors=FALSE)
  
  #Single copy
  orthogroups.stats$copy_number[match(rownames(orthogroups.copies.ingroup[apply(orthogroups.copies.ingroup < 2, 1, all),]), orthogroups.stats$orthogroup)] <- "single"
  
  #Secretome type
  
  print("Assigning orthogroups as core, accessory or specific")
  #Create bar to show progress
  progress.bar <- txtProgressBar(1, length(rownames(orthogroups.copies.ingroup)), initial=0, char="=", style=3)
  for (j in 1:length(rownames(orthogroups.copies.ingroup))) {
    
    #Update progress bar
    setTxtProgressBar(progress.bar, j)
    
    if (length(which(orthogroups.copies.ingroup[j,] == 0)) == (length(colnames(orthogroups.copies.ingroup)) - 1)) {
      orthogroups.stats$secretome[orthogroups.stats$orthogroup == rownames(orthogroups.copies.ingroup)[j]] <- "specific"
    }
    if (length(which(orthogroups.copies.ingroup[j,] == 0)) == 0) {
      orthogroups.stats$secretome[orthogroups.stats$orthogroup == rownames(orthogroups.copies.ingroup)[j]] <- "core"
    }
    if (length(which(orthogroups.copies.ingroup[j,] == 0)) < (length(colnames(orthogroups.copies.ingroup)) - 1) && length(which(orthogroups.copies.ingroup[j,] == 0)) > 0) {
      orthogroups.stats$secretome[orthogroups.stats$orthogroup == rownames(orthogroups.copies.ingroup)[j]] <- "accessory"
    }
  }
  close(progress.bar)
  
  #Effectors
  
  #Assess which orthogroups contain only effectors (SP-only) or both effectors and other genes (SP-mixed) 
  effectors <- names(which(rowSums(effector.count.ingroup) > 0))
  
  effector.type <- list()
  
  #For each sample...
  for (i in 1:length(colnames(orthogroups.copies.ingroup))) {
    
    #Print progress
    cat("Assigning effector orthogroups as effector-only or effector-mixed: ", (i - 1), "/", length(colnames(orthogroups.copies.ingroup)), " taxa", "\r")
    
    for (j in effectors) {
      if (orthogroups.copies.ingroup[j,i] == effector.count.ingroup[j,i]) {
        effector.type[[j]][i] <- "effector"
      } else {
        effector.type[[j]][i] <- "mixed"
      }
    }
    
  }
  print(paste0("Assigning effector orthogroups as effector-only or effector-mixed: ", i, "/", length(colnames(orthogroups.copies.ingroup)), " taxa"))
  
  for (j in effectors) {
    if(grepl("mixed", effector.type[j]) == TRUE) {
      orthogroups.stats$effector[orthogroups.stats$orthogroup == j] <- "mixed"
    } else {
      orthogroups.stats$effector[orthogroups.stats$orthogroup == j] <- "effector"
    }
  }

  assign(paste0("orthogroups.stats.ingroup", ingroup), orthogroups.stats)
  assign(paste0("orthogroups.copies.ingroup", ingroup), orthogroups.copies.ingroup)
  assign(paste0("effector.count.ingroup", ingroup), effector.count.ingroup)
  
  #print(paste0("Summary of orthogroups saved in orthogroups_stats-", Sys.Date(), ".csv"))
  #write.csv(orthogroups.stats,
  #          file=paste0("orthogroups_stats-", Sys.Date(), ".csv"), row.names=FALSE, quote=FALSE)
  
}

#Lists of orthogroups for selection analyses
#Core, single-copy effectors (for aBSREL and MEME)
core.SC.effectors <- Reduce(intersect,
                            list(orthogroups.stats.ingroup0$orthogroup[which(orthogroups.stats.ingroup0$copy_number == "single")],
                                 orthogroups.stats.ingroup0$orthogroup[which(orthogroups.stats.ingroup0$secretome == "core")],
                                 orthogroups.stats.ingroup0$orthogroup[which(orthogroups.stats.ingroup0$effector == "effector")]))
#Core and accessory (in >=3 taxa), single-copy orthogroups (for BUSTED)
accessory.core.SC <- Reduce(intersect, list(names(which(rowSums(orthogroups.copies > 0) >= 3)),
                                            orthogroups.stats.ingroup1$orthogroup[which(orthogroups.stats.ingroup1$copy_number == "single")]))

write.table(accessory.core.SC,
            file="../selection/orthogroups_busted.csv", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(core.SC.effectors,
            file="../selection/orthogroups_meme.csv", col.names=FALSE, row.names=FALSE, quote=FALSE)


print("Orthogroups for selection analyses saved in selection directory")

print(paste0("Workspace saved in effector-matrices-", Sys.Date(), ".RData"))
save.image(file=paste0("effector-matrices-", Sys.Date(), ".RData"))
