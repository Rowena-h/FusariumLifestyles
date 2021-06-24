#Script to estimate scale parameter (beta) for the MCMCTree rate prior

library(ape)
library(adephylo)

#Read in and root species tree
tree <- read.tree("../../species_tree/iqtree/gene_partitions/fus_proteins_62T_iqtree_genepart.contree")
tree <- root(tree, outgroup="Ilysp1_GeneCatalog_proteins_20121116", edgelabel=TRUE, resolve.root=TRUE)

#beta = (alpha x root-time) / mean tip to root distance
1 / mean(distRoot(tree, tree$tip.label, method="patristic"))
