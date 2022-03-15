#!/usr/bin/env Rscript
#Script to estimate scale parameter (beta) for the MCMCTree rate prior

library(ape)
library(adephylo)

#Read in rooted species trees
raxmlng <- read.tree("../fus_proteins_62T.raxml.support_rooted")
iqtree <- read.tree("../../phylogenomics/species_tree/iqtree/fus_proteins_62T_iqtree_genepart.contree")
iqtree <- root(iqtree, outgroup="Ilysp1_GeneCatalog_proteins_20121116", edgelabel=TRUE, resolve.root=TRUE)

#Calculate substitution rate prior
#beta = (alpha x root-time) / mean tip to root distance
print(mean(distRoot(iqtree, iqtree$tip.label, method="patristic")))
print(mean(distRoot(raxmlng, raxmlng$tip.label, method="patristic")))
print(paste0("IQ-TREE: ", 1 / mean(distRoot(iqtree, iqtree$tip.label, method="patristic"))))
print(paste0("RAxML-NG: ", 1 / mean(distRoot(raxmlng, raxmlng$tip.label, method="patristic"))))
