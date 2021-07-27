#!/usr/bin/env Rscript
##Script to prepare trees for Contrast-FEL##

library(ape)

#Read in metadata
metadata <- read.csv("../metadata.csv")

#List of SC core gene trees
trees <- list.files("trees/", pattern=".bestTree_rooted$")

#For each tree...
for (i in trees) {
  
  #Read in tree
  tree <- read.tree(paste0("trees/", i))
  
  #For each lifestyle...
  for (j in unique(metadata$lifestyle)) {
    
    #Add label to taxa in that lifestyle
    tree.edit <- tree
    tree.edit$tip.label[match(metadata$tip[metadata$lifestyle == j], tree.edit$tip.label)] <- paste0(tree.edit$tip.label[match(metadata$tip[metadata$lifestyle == j], tree.edit$tip.label)], "{lifestyle}")
    
    #Write tree
    lifestyle <- sub(" ", "", j)
    write.tree(tree.edit, paste0("trees/", i, ".", lifestyle))
    
  }
  
}
