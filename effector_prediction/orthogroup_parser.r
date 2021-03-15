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

#Assess which orthogroups contain only SPs (SP-only) or both SPs and other genes (SP-mixed) 
mixed <- list()

#For each sample...
for (i in 1:length(colnames(orthogroups.copies))) {
  
  #Print progress
  cat("Assigning orthogroups as SP-only or SP-mixed: ", (i - 1), "/", length(colnames(orthogroups.copies)), " taxa", "\r")
  
  for (j in 1:length(rownames(orthogroups.copies))) {
    if (orthogroups.copies[j,i] == effector.count[j,i]) {
      mixed[[rownames(orthogroups.copies)[j]]][i] <- "SP-only"
    } else {
      mixed[[rownames(orthogroups.copies)[j]]][i] <- "SP-mixed"
    }
  }
  
}
print(paste0("Assigning orthogroups as SP-only or SP-mixed: ", i, "/", length(colnames(orthogroups.copies)), " taxa"))

#Assess which orthogroups are in all species (core), one species (specific), or some species (accessory)
secretome <- vector(mode="character", length=length(rownames(effector.count)))

print("Assigning orthogroups as core, accessory or specific")
#Create bar to show progress
progress.bar <- txtProgressBar(1, length(rownames(effector.count)), initial=0, char="=", style=3)
for (j in 1:length(rownames(effector.count))) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, j)
  
  if (length(which(effector.count[j,] == 0)) == (length(colnames(effector.count)) - 1)) {
    secretome[j] <- "specific"
  }
  if (length(which(effector.count[j,] == 0)) == 0) {
    secretome[j] <- "core"
  }
  if (length(which(effector.count[j,] == 0)) < (length(colnames(effector.count)) - 1) && length(which(effector.count[j,] == 0)) > 0) {
    secretome[j] <- "accessory"
  }
}
close(progress.bar)

#Add column to effector count dataframe with whether orthogroup is SP-only or SP-mixed
effector.count$mixed <- NA

for (j in 1:length(rownames(effector.count))) {
  if(grepl("SP-mixed", mixed[j]) == TRUE) {
    effector.count$mixed[j] <- "SP-mixed"
  } else {
    effector.count$mixed[j] <- "SP-only"
  }
}

#Add column to effector count dataframe with whether orthogroup core, specific or accessory
effector.count$secretome <- secretome
#Reorder dataframe columns
effector.count <- effector.count %>% select(secretome, mixed, everything())

#Filter out non SP orthogroups
effector.count <- effector.count[effector.count$secretome != "",]
#Filter SP-only orthogroups
effector.count.SP <- effector.count[effector.count$mixed == "SP-only",]
#Filter for single-copy orthogroups
effector.count.SC <- effector.count[apply(effector.count[3:length(colnames(effector.count))] < 2, 1, all),]
#Filter for single-copy, SP-only orthogroups
effector.count.SC.SP <- effector.count.SC[effector.count.SC$mixed == "SP-only",]


#List of core SP-only single-copy orthogroups to check for positive selection
selection <- rownames(effector.count.SC.SP[effector.count.SC.SP$secretome == "core",])
#List of core SP-only multi-copy orthogroups to check for positive selection
selection.multi <- rownames(effector.count.SP[effector.count.SP$secretome == "core",])
selection.multi <- setdiff(selection.multi, selection)
         
write.table(selection,
            file="selection/orthogroups_selection.csv", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(selection.multi,
            file="selection/orthogroups_selection_multi.csv", col.names=FALSE, row.names=FALSE, quote=FALSE)


print("Core effector orthogroups saved in selection/orthogroups_selection.csv")