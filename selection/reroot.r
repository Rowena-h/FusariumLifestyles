#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

#Test if there are two arguments: if not, return an error
if (length(args)<2) {
  stop("Two arguments must be supplied: a tree first and an outgroup string second", call.=FALSE)
} 

library(ape)

tree <- read.tree(args[1])
outgroup <- as.character(args[2])

tree <- root(tree, outgroup, resolve.root=TRUE, edgelabel=TRUE)
write.tree(tree, file=paste0(args[1], "_rooted"))


