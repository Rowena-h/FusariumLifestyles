#!/usr/bin/env Rscript
#Script to estimate scale parameter (beta) for the MCMCTree rate prior

library(ape)
library(adephylo)

#Read in rooted species trees
iqtree <- read.tree("fus_proteins_62T_iqtree_genepart.contree_rooted")
raxmlng <- read.tree("../../phylogenomics/species_tree/raxml-ng/fus_proteins_62T.raxml.support")
raxmlng <- root(raxmlng, outgroup="Ilysp1_GeneCatalog_proteins_20121116", edgelabel=TRUE, resolve.root=TRUE)

#Calculate substitution rate prior
#beta = (alpha x root-time) / mean tip to root distance
print(paste0("IQ-TREE: ", 1 / mean(distRoot(iqtree, tree$tip.label, method="patristic"))))
print(paste0("RAxML-NG: ", 1 / mean(distRoot(raxmlng, tree$tip.label, method="patristic"))))
