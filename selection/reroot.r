#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

#Test if there are two arguments: if not, return an error
if (length(args)<3) {
  stop("Three arguments must be supplied: a tree first, an outgroup string second and the directory to save the rooted tree third.", call.=FALSE)
} 

library(ape)

tree <- read.tree(args[1])
outgroup <- as.character(args[2])

file <- sub(".*\\/", "", args[1])

tree <- root(tree, outgroup, resolve.root=TRUE, edgelabel=TRUE)
write.tree(tree, file=paste0(args[3], file, "_rooted"))


