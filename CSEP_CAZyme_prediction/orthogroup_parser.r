#!/usr/bin/env Rscript
##Script to parse OrthoFinder, CSEPfilter and run_dbcan results into presence-absence matrix of orthogroups for each taxon##

library(dplyr)
library(rvest)

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
message("Reading in CSEPs and CAZymes")
for (i in colnames(orthogroups)) {
  
  #Read in CSEP tables
  CSEPs <- read.table(paste0(i, ".faa_cseps"), sep="\t", header=FALSE,
                      col.names=c("Gene_ID", "Protein_ID", "PHI-base_entry",
                                  "Gene_name", "Pathogen_ID", "Organism", "Mutant_phenotype",
                                  "e-value", "bitscore", "pident"), fill=TRUE,
                      comment.char="")
  
  #Read in CAZyme tables
  CAZymes <- read.table(paste0(i, ".faa_cazymes"), sep="\t", header=FALSE,
                        col.names=c("Gene_ID", "EC#", "HMMER",
                                    "eCAMI", "DIAMOND", "#ofTools"),
                        comment.char="")
  
  #Replace pipes (|) with hyphens
  CSEPs$Gene_ID <- gsub("\\|", "-", CSEPs$Gene_ID)
  CAZymes$Gene_ID <- gsub("\\|", "-", CAZymes$Gene_ID)
  orthogroups[,i] <- gsub("\\|", "-", orthogroups[,i])
  
  assign(paste0(i, ".CSEPs"), CSEPs)
  assign(paste0(i, ".CAZymes"), CAZymes)
  
}

#Make dataframe with gene counts for each orthogroup
orthogroups.count <- orthogroups

for (i in 1:length(colnames(orthogroups.count))) {
  orthogroups.count[, i] <- sapply(strsplit(orthogroups[, i], " "), length)
}

#Make dataframe for CSEP/CAZyme counts
count <- data.frame(matrix(0, ncol=ncol(orthogroups), nrow=nrow(orthogroups)))
colnames(count) <- colnames(orthogroups)
rownames(count) <- rownames(orthogroups)


message("Counting number of CSEPs/CAZymes in each orthogroup:")

#For CSEPs and CAZymes in turn...
for (group in c("CSEP", "CAZyme")) {
  
  message(group, "s:")
  tmp.count <- count
  #Make list to capture metadata
  gene.list <- list()
  
  #For each sample...
  for (i in 1:length(colnames(tmp.count))) {
    
    #Print progress
    message((i - 1), "/", length(colnames(tmp.count)))
    
    #Retrieve the list of gene predictions
    genes <- get(paste0(colnames(tmp.count)[i], ".", group, "s"))
    
    #For each gene...
    for (j in 1:length(genes$Gene_ID)) {
      
      #Retrieve gene
      gene <- grep(genes$Gene_ID[j], orthogroups[, i])
      
      #Search for gene in corresponding column of orthogroups dataframe and add 1 to orthogroup count
      tmp.count[gene, i] <- tmp.count[gene, i] + 1
      
      #Get metadata
      gene.list[[rownames(tmp.count[gene,])]][[colnames(tmp.count)[i]]] <- genes[j,]
      
    }
    
  }
  
  message(i, "/", length(colnames(tmp.count)))
  assign(paste0(group, ".count"), tmp.count)
  assign(paste0(group, ".list"), gene.list)
  
}


#Read in sample metadata
metadata <- read.csv("../metadata.csv")

#Including and excluding the outgroup...
for (ingroup in c(0, 1)) {
  
  if (ingroup == 1) {
    
    #Filter for taxa in ingroup
    orthogroups.count.ingroup <- orthogroups.count[which(!is.na(match(colnames(orthogroups.count),
                                                                      metadata$file2[metadata$ingroup == ingroup])))]
    CSEP.count.ingroup <- CSEP.count[which(!is.na(match(colnames(CSEP.count),
                                                        metadata$file2[metadata$ingroup == ingroup])))]
    CAZyme.count.ingroup <- CAZyme.count[which(!is.na(match(colnames(CAZyme.count),
                                                            metadata$file2[metadata$ingroup == ingroup])))]
    
  } else {
    
    orthogroups.count.ingroup <- orthogroups.count
    CSEP.count.ingroup <- CSEP.count
    CAZyme.count.ingroup <- CAZyme.count
    
  }
  
  #Make dataframe to collect stats for orthogroups
  orthogroups.stats <- data.frame(orthogroup=rownames(orthogroups.count.ingroup),
                                  copy_number="multi",
                                  category=NA,
                                  CSEP=NA,
                                  CAZyme=NA,
                                  stringsAsFactors=FALSE)
  
  #Add whether orthogroups are single copy
  orthogroups.stats$copy_number[match(rownames(orthogroups.count.ingroup[apply(orthogroups.count.ingroup < 2, 1, all),]), orthogroups.stats$orthogroup)] <- "single"
  
  #Categorise orthogroups according to sharedness
  message("Assigning orthogroups as core, accessory or specific for ingroup ", ingroup, ":")
  #Create bar to show progress
  progress.bar <- txtProgressBar(1, length(rownames(orthogroups.count.ingroup)), initial=0, char="=", style=3)
  for (j in 1:length(rownames(orthogroups.count.ingroup))) {
    
    #Update progress bar
    setTxtProgressBar(progress.bar, j)
    
    if (length(which(orthogroups.count.ingroup[j,] == 0)) == (length(colnames(orthogroups.count.ingroup)) - 1)) {
      orthogroups.stats$category[orthogroups.stats$orthogroup == rownames(orthogroups.count.ingroup)[j]] <-
        "specific"
    }
    if (length(which(orthogroups.count.ingroup[j,] == 0)) == 0) {
      orthogroups.stats$category[orthogroups.stats$orthogroup == rownames(orthogroups.count.ingroup)[j]] <-
        "core"
    }
    if (length(which(orthogroups.count.ingroup[j,] == 0)) < (length(colnames(orthogroups.count.ingroup)) - 1) && length(which(orthogroups.count.ingroup[j,] == 0)) > 0) {
      orthogroups.stats$category[orthogroups.stats$orthogroup == rownames(orthogroups.count.ingroup)[j]] <-
        "accessory"
    }
  }
  close(progress.bar)
  
  #Determine which orthogroups contain only CSEPs (SP-only) or both CSEPs and other genes (SP-mixed) 
  CSEPs <- names(which(rowSums(CSEP.count.ingroup) > 0))
  CAZymes <- names(which(rowSums(CAZyme.count.ingroup) > 0))
  
  CSEP.type <- list()
  CAZyme.type <- list()
  
  message("Assigning CSEP/CAZyme orthogroups as -only or -mixed for ingroup ", ingroup, ":")
  
  #For each sample...
  for (i in 1:length(colnames(orthogroups.count.ingroup))) {
    
    #Print progress
    message((i - 1), "/", length(colnames(orthogroups.count.ingroup)))
    
    for (j in CSEPs) {
      if (orthogroups.count.ingroup[j,i] == CSEP.count.ingroup[j,i]) {
        CSEP.type[[j]][i] <- "CSEP"
      } else {
        CSEP.type[[j]][i] <- "mixed"
      }
    }
    
    for (k in CAZymes) {
      if (orthogroups.count.ingroup[k,i] == CAZyme.count.ingroup[k,i]) {
        CAZyme.type[[k]][i] <- "CAZyme"
      } else {
        CAZyme.type[[k]][i] <- "mixed"
      }
    }
    
  }
  message(i, "/", length(colnames(orthogroups.count.ingroup)))
  
  for (j in CSEPs) {
    if(grepl("mixed", CSEP.type[j]) == TRUE) {
      orthogroups.stats$CSEP[orthogroups.stats$orthogroup == j] <- "mixed"
    } else {
      orthogroups.stats$CSEP[orthogroups.stats$orthogroup == j] <- "CSEP"
    }
  }
  
  for (j in CAZymes) {
    if(grepl("mixed", CAZyme.type[j]) == TRUE) {
      orthogroups.stats$CAZyme[orthogroups.stats$orthogroup == j] <- "mixed"
    } else {
      orthogroups.stats$CAZyme[orthogroups.stats$orthogroup == j] <- "CAZyme"
    }
  }
  
  assign(paste0("orthogroups.stats.ingroup", ingroup), orthogroups.stats)
  assign(paste0("orthogroups.count.ingroup", ingroup), orthogroups.count.ingroup)
  assign(paste0("CSEP.count.ingroup", ingroup), CSEP.count.ingroup)
  assign(paste0("CAZyme.count.ingroup", ingroup), CAZyme.count.ingroup)
  
}


message("Adding CSEP/CAZyme annotations:")

orthogroups.stats <- data.frame(orthogroups.stats,
                                CSEP_name=NA, PHI.base_entry=NA,
                                EC=NA, CAZyme_family=NA, CAZyme_name=NA,
                                PCWDE=NA)

#Create bar to show progress
progress.bar <- txtProgressBar(1, length(names(CAZyme.list)), initial=0, char="=", style=3)

#For each CAZyme...
for (i in 1:length(names(CAZyme.list))) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, i)
  
  #Get EC numbers
  EC.number <- unique(unlist(lapply(CAZyme.list[[i]], "[", "EC.")))
  
  #If gene is classified as more than one EC number...
  if (length(EC.number > 1)) {
    
    #Search for each EC number separately
    for (j in EC.number) {
      
      j <- sub("\\|.*", "", j)
      
      enzyme.list <- list()
      
      #Search ExplorEnz website for EC number
      explorenz <- read_html(paste0("https://www.enzyme-database.org/query.php?ec=", j))
      #Scrape accepted name
      enzyme.list[EC.number] <- explorenz %>%
        html_element("tr:nth-child(2) td+ td") %>% 
        html_text2()
      
    }
    
    #Summarise EC numbers/accepted names
    EC.number <- paste(EC.number, collapse=",")
    CAZyme.name <- paste(unique(unlist(enzyme.list)), collapse=",")
    
  } else {
    
    #Search ExplorEnz website for EC number
    explorenz <- read_html(paste0("https://www.enzyme-database.org/query.php?ec=", EC.number))
    #Scrape accepted name
    CAZyme.name <- explorenz %>%
      html_element("tr:nth-child(2) td+ td") %>% 
      html_text2()
    
  }
  
  #Add results to dataframe
  orthogroups.stats$EC[match(names(CAZyme.list)[i], orthogroups.stats$orthogroup)] <- EC.number
  orthogroups.stats$CAZyme_name[match(names(CAZyme.list)[i], orthogroups.stats$orthogroup)] <- CAZyme.name
  
  #Get CAZyme family
  CAZyme.family <- unique(unlist(lapply(CAZyme.list[[i]], "[", "DIAMOND")))
  
  #If gene is classified as more than one CAZyme family...
  if (length(CAZyme.family > 1)) {
    
    #Summarise families
    CAZyme.family <- paste(CAZyme.family, collapse=",")
    
  }
  
  #Add to dataframe
  orthogroups.stats$CAZyme_family[match(names(CAZyme.list)[i], orthogroups.stats$orthogroup)] <- CAZyme.family
  
}
close(progress.bar)


#Create bar to show progress
progress.bar <- txtProgressBar(1, length(names(CSEP.list)), initial=0, char="=", style=3)

#For each CSEP...
for (i in 1:length(names(CSEP.list))) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, i)
  
  #Get gene name
  CSEP.name <- unique(unlist(lapply(CSEP.list[[i]], "[", "Gene_name")))
  #Get PHI-base entry
  PHI.base <- unique(unlist(lapply(CSEP.list[[i]], "[", "PHI.base_entry")))
  
  #If there's an associated name...
  if (CSEP.name != "") {
    
    #If gene is classified under more than one gene name...
    if (length(CSEP.name > 1)) {
      
      #Summarise names
      CSEP.name <- paste(CSEP.name[CSEP.name != ""], collapse=",")
      PHI.base <- paste(PHI.base[PHI.base != ""], collapse=",")
      
    }
    
    #Add results to dataframe
    orthogroups.stats$CSEP_name[match(names(CSEP.list)[i], orthogroups.stats$orthogroup)] <- CSEP.name
    orthogroups.stats$PHI.base_entry[match(names(CSEP.list)[i], orthogroups.stats$orthogroup)] <- PHI.base
    
  }
  
}
close(progress.bar)

#Add to all dataframes
orthogroups.stats.ingroup0 <- 
  data.frame(orthogroups.stats.ingroup0,
             orthogroups.stats[match(orthogroups.stats.ingroup0$orthogroup,
                                     orthogroups.stats$orthogroup),
                               c("CSEP_name", "PHI.base_entry", "EC",
                                 "CAZyme_family", "CAZyme_name", "PCWDE")])

orthogroups.stats.ingroup1 <- 
  data.frame(orthogroups.stats.ingroup1,
             orthogroups.stats[match(orthogroups.stats.ingroup1$orthogroup,
                                     orthogroups.stats$orthogroup),
                               c("CSEP_name", "PHI.base_entry", "EC",
                                 "CAZyme_family", "CAZyme_name", "PCWDE")])

#Remove orthogroups missing from ingroup
orthogroups.stats.ingroup1 <- orthogroups.stats.ingroup1[!is.na(orthogroups.stats.ingroup1$category),]

message(paste0("Results saved in orthogroup-matrices-", Sys.Date(), ".RData"))
save(orthogroups.stats.ingroup0,
     orthogroups.stats.ingroup1,
     orthogroups.count.ingroup0,
     orthogroups.count.ingroup1,
     CSEP.count.ingroup0,
     CSEP.count.ingroup1,
     CAZyme.count.ingroup0,
     CAZyme.count.ingroup1,
     file=paste0("orthogroup-matrices-", Sys.Date(), ".RData"))
