#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

#Test if there are two arguments: if not, return an error
if (length(args)<2) {
  stop("Two arguments must be supplied: a tree first and the directory to save the rooted tree second.", call.=FALSE)
} 

library(ape)

#Read in tree
tree <- read.tree(args[1])

file <- sub(".*\\/", "", args[1])

#Remove branch lengths and node and tip labels
tree$edge.length <- NULL
tree$node.label <- NULL
tree$tip.label <- substr(tree$tip.label, 1, 13) 

write.tree(tree, file=paste0(args[2], file, "_blank"))


