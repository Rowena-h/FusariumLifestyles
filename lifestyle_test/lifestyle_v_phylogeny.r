##https://github.com/fantin-mesny/Effect-Of-Biological-Categories-On-Genomes-Composition

library(ape)
library(dplyr)

load("../effector_prediction/effector-matrices-2021-05-25.RData")

##LIFESTYLE VERSUS PHYLOGENY TEST

metadata <- read.csv("../metadata.csv")

#Read in tree
phy <- read.tree("../phylogenomics/species_tree/iqtree/gene_partitions/fus_proteins_62T_iqtree_genepart.contree")
#Root tree
outgroup <- "Ilysp1_GeneCatalog_proteins_20121116"
phy <- root(phy, outgroup, resolve.root=TRUE, edgelabel=TRUE)
#Remove outgroup from tree and write to file
phy.ingroup <- drop.tip(phy, outgroup)
write.tree(phy.ingroup, "species_tree_ingroup.tre")

#Effectors
#Transpose count dataframe (excluding outgroup)
lifestyle.test.effectors <- as.data.frame(t(effector.count[-which(rowSums(effector.count) == 0),][-which(colnames(effector.count) == outgroup)]))
#Add column with names
lifestyle.test.effectors$genome <- rownames(lifestyle.test.effectors)
#Add column with lifestyle
lifestyle.test.effectors$lifestyle <- gsub(" ", "", metadata$lifestyle[match(rownames(lifestyle.test.effectors), metadata$file)])
#Replace spaces and hyphens to match tree labels
lifestyle.test.effectors$lifestyle <- gsub("-", "", lifestyle.test.effectors$lifestyle)
lifestyle.test.effectors$genome <- gsub("\\.", "_", lifestyle.test.effectors$genome)
lifestyle.test.effectors$genome <- gsub(" ", "_", lifestyle.test.effectors$genome)
#Reorder columns
lifestyle.test.effectors <- lifestyle.test.effectors %>% select(genome, lifestyle, everything())

#Write to file
write.csv(lifestyle.test.effectors, "lifestyle-test-effectors.csv", row.names=FALSE, quote=FALSE)

#Repeat for orthogroups
lifestyle.test.orthogroups <- as.data.frame(t(orthogroups.copies[-which(colnames(orthogroups.copies) == outgroup)]))
lifestyle.test.orthogroups$genome <- rownames(lifestyle.test.orthogroups)
lifestyle.test.orthogroups$lifestyle <- gsub(" ", "", metadata$lifestyle[match(rownames(lifestyle.test.orthogroups), metadata$file)])
lifestyle.test.orthogroups$lifestyle <- gsub("-", "", lifestyle.test.orthogroups$lifestyle)
lifestyle.test.orthogroups$genome <- gsub("\\.", "_", lifestyle.test.orthogroups$genome)
lifestyle.test.orthogroups$genome <- gsub(" ", "_", lifestyle.test.orthogroups$genome)
lifestyle.test.orthogroups <- lifestyle.test.orthogroups %>% select(genome, lifestyle, everything())

write.csv(lifestyle.test.orthogroups, "lifestyle-test-orthogroups.csv", row.names=FALSE, quote=FALSE)

system("qsub lifestyle-test.sh")
