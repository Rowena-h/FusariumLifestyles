#Script to estimate scale parameter (beta) for the MCMCTree rate prior

library(ape)
library(adephylo)

#Read in and root ML species trees
iqtree <- read.tree("../../phylogenomics/species_tree/iqtree/fus_proteins_62T_iqtree_genepart.contree")
iqtree <- root(iqtree, outgroup="Ilysp1_GeneCatalog_proteins_20121116", edgelabel=TRUE, resolve.root=TRUE)

raxmlng <- read.tree("../../phylogenomics/species_tree/raxml-ng/fus_proteins_62T.raxml.support")
raxmlng <- root(raxmlng, outgroup="Ilysp1_GeneCatalog_proteins_20121116", edgelabel=TRUE, resolve.root=TRUE)

#beta = (alpha x root-time) / mean tip to root distance
print(paste0("IQ-TREE: ", 1 / mean(distRoot(iqtree, tree$tip.label, method="patristic"))))
print(paste0("RAxML-NG: ", 1 / mean(distRoot(raxmlng, tree$tip.label, method="patristic"))))

