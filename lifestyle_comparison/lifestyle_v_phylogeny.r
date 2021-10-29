#!/usr/bin/env Rscript
##Script to produce input files for lifestyle test:
##https://github.com/fantin-mesny/Effect-Of-Biological-Categories-On-Genomes-Composition

library(ape)
library(MCMCtreeR)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

#Test if there is one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied: the orthology parsing results RData file to be used", call.=FALSE)
} 

#Load orthogroup presence absence matrices
load(args[1])

##LIFESTYLE VERSUS PHYLOGENY TEST

metadata <- read.csv("../metadata.csv")

#Read in dated tree
phy <- readMCMCtree("../divergence_time_estimation/mcmctree/run1_independent/FigTree.tre", 
		    forceUltrametric=TRUE)$apePhy
#Remove outgroup from tree and write to file
outgroup <- "Ilysp1_GeneCa"
phy.ingroup <- drop.tip(phy, outgroup)
write.tree(phy.ingroup, "species_tree_ingroup.tre")

#CSEPs
#Transpose count dataframe (excluding outgroup)
lifestyle.test.CSEPs <- as.data.frame(t(CSEP.count.ingroup1))
#Add column with names
lifestyle.test.CSEPs$genome <- rownames(lifestyle.test.CSEPs)
#Add column with lifestyle
lifestyle.test.CSEPs$lifestyle <- gsub(" ", "", metadata$lifestyle[match(rownames(lifestyle.test.CSEPs), metadata$file2)])
#Replace spaces and hyphens to match tree labels
lifestyle.test.CSEPs$genome <- substr(lifestyle.test.CSEPs$genome, 1, 13)
#Reorder columns
lifestyle.test.CSEPs <- lifestyle.test.CSEPs %>% select(genome, lifestyle, everything())

#Write to file
write.csv(lifestyle.test.CSEPs, "lifestyle-test-CSEPs.csv", row.names=FALSE, quote=FALSE)

#Repeat for orthogroups
lifestyle.test.orthogroups <- as.data.frame(t(orthogroups.copies.ingroup1))
lifestyle.test.orthogroups$genome <- rownames(lifestyle.test.orthogroups)
lifestyle.test.orthogroups$lifestyle <- gsub(" ", "", metadata$lifestyle[match(rownames(lifestyle.test.orthogroups), metadata$file2)])
lifestyle.test.orthogroups$genome <- substr(lifestyle.test.orthogroups$genome, 1, 13)
lifestyle.test.orthogroups <- lifestyle.test.orthogroups %>% select(genome, lifestyle, everything())

write.csv(lifestyle.test.orthogroups, "lifestyle-test-orthogroups.csv", row.names=FALSE, quote=FALSE)

#Submit test
system("qsub lifestyle-test.sh")
