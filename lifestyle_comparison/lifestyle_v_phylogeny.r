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

for (i in c("orthogroups", "CSEP", "CAZyme")) {
  
  #Transpose count dataframe (excluding outgroup)
  lifestyle.test <- as.data.frame(t(get(paste0(i, ".count.ingroup1"))))
  #Add column with names
  lifestyle.test$genome <- rownames(lifestyle.test)
  #Add column with lifestyle
  lifestyle.test$lifestyle <- gsub(" ", "", metadata$lifestyle[match(rownames(lifestyle.test), metadata$file2)])
  #Replace names to match dated tree labels
  lifestyle.test$genome <- metadata$short.tip[match(lifestyle.test$genome, metadata$file2)]
  #Reorder columns
  lifestyle.test <- lifestyle.test %>% select(genome, lifestyle, everything())
  
  #Write to file
  write.csv(lifestyle.test, paste0("lifestyle-test-", i, ".csv"), row.names=FALSE, quote=FALSE)
  
}

#Submit test
system("qsub lifestyle-test.sh")
