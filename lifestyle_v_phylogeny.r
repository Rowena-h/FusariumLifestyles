setwd("D:/Documents/GitHub/FusariumEndophytes/")
##https://github.com/fantin-mesny/Effect-Of-Biological-Categories-On-Genomes-Composition

library(ape)

##LIFESTYLE VERSUS PHYLOGENY TEST

metadata <- read.csv("metadata.csv")

#Read in tree
phy <- read.tree("phylogenomics/species_tree/astral/fus_astral_proteins_62T.tre")
#Root tree
outgroup <- "Ilysp1_GeneCatalog_proteins_20121116"
phy <- root(phy, outgroup, resolve.root=TRUE, edgelabel=TRUE)
#Remove outgroup from tree and write to file
phy.ingroup <- drop.tip(phy, outgroup)
write.tree(phy.ingroup, "phylogenomics/species_tree/RAxML_bipartitions.fus_proteins_27T_bsadd_rooted_ingroup")

#Replace sample IDs with names
#colnames(effector.count.SC.SP) <- c("secretome", "mixed", metadata$name[match(colnames(effector.count.SC.SP)[-c(1, 2)], metadata$file)])
#colnames(orthogroups.copies) <- metadata$name[match(colnames(orthogroups.copies), metadata$file)]

#Effectors
#Transpose count dataframe (excluding outgroup)
lifestyle.test.effectors <- as.data.frame(t(effector.count.SC.SP[-c(1, 2, which(colnames(effector.count.SC.SP) == outgroup))]))
#Add column with names
lifestyle.test.effectors$genome <- rownames(lifestyle.test.effectors)
#Add column with lifestyle
lifestyle.test.effectors$lifestyle <- gsub(" ", "", metadata$lifestyle.hyp1[match(rownames(lifestyle.test.effectors), metadata$file)])
#Replace spaces and hyphens to match tree labels
lifestyle.test.effectors$lifestyle <- gsub("-", "", lifestyle.test.effectors$lifestyle)
lifestyle.test.effectors$genome <- gsub("\\.", "_", lifestyle.test.effectors$genome)
lifestyle.test.effectors$genome <- gsub(" ", "_", lifestyle.test.effectors$genome)
#Reorder columns
lifestyle.test.effectors <- lifestyle.test.effectors[c(length(colnames(lifestyle.test.effectors))-1,
                                                       length(colnames(lifestyle.test.effectors)),
                                                       1:(length(colnames(lifestyle.test.effectors))-2))]
#Write to file
#write.csv(lifestyle.test.effectors, "lifestyle-test-effectors.csv", row.names=FALSE)

#Repeat for orthogroups
lifestyle.test.orthogroups <- as.data.frame(t(orthogroups.copies[-which(colnames(orthogroups.copies) == outgroup)]))
lifestyle.test.orthogroups$genome <- rownames(lifestyle.test.orthogroups)
lifestyle.test.orthogroups$lifestyle <- gsub(" ", "", metadata$lifestyle.hyp1[match(rownames(lifestyle.test.orthogroups), metadata$file)])
lifestyle.test.orthogroups$lifestyle <- gsub("-", "", lifestyle.test.orthogroups$lifestyle)
lifestyle.test.orthogroups$genome <- gsub("\\.", "_", lifestyle.test.orthogroups$genome)
lifestyle.test.orthogroups$genome <- gsub(" ", "_", lifestyle.test.orthogroups$genome)
lifestyle.test.orthogroups <- lifestyle.test.orthogroups[c(length(colnames(lifestyle.test.orthogroups))-1,
                                                           length(colnames(lifestyle.test.orthogroups)),
                                                           1:(length(colnames(lifestyle.test.orthogroups))-2))]
#write.csv(lifestyle.test.orthogroups, "lifestyle-test-orthogroups.csv", row.names=FALSE)