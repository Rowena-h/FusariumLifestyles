setwd("D:/Documents/Kew/Ascomycota/Musa Fusaria/Comparative genomics/PredGPI/")

#Retrieve GPI anchor predictions from PredGPI (http://gpcr.biocomp.unibo.it/predgpi/)

library(ragp)
library(seqinr)

list <- read.csv("fus_list", header=FALSE)
list <- list$V1
  
for (i in list) {
  #Read in protein sets for each assembly
  tmp <- read.fasta(i, seqtype="AA", as.string=TRUE)
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
