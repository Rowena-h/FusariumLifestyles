#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

#Test if there are two arguments: if not, return an error
if (length(args)<2) {
  stop("Two arguments must be supplied: the directory where the protein fasta files are located first and a file listing the sample names second", call.=FALSE)
} 

#Retrieve GPI anchor predictions from PredGPI (http://gpcr.biocomp.unibo.it/predgpi/)

library(ragp)
library(seqinr)

list <- read.csv(paste0("../", args[2]), header=FALSE)
list <- list$V1
  
for (i in list) {
  #Read in protein sets for each assembly
  tmp <- read.fasta(paste0(args[1], i), seqtype="AA", as.string=TRUE)
  #Replace protein IDs as limited to 30 characters
  df <- data.frame(orig=names(tmp), replace=seq(length(names(tmp))))
  names(tmp) <- df$replace
  results <- get_pred_gpi(tmp, spec=0.99, progress=TRUE)
  #Return original protein IDs
  results$id <- df$orig
  #Print list of GPI-anchored proteins
  gpilist <- df$orig[results$is.gpi == TRUE]
  assign(paste0(i,"_gpilist"), gpilist)
  #Write file
  write(gpilist, file=paste0(i, "_predgpi_GPlist"))
}
