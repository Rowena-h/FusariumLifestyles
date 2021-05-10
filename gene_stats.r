library(seqinr)
library(ggplot2)
library(jsonlite)

taxa <- list.files("orthology_inference/", pattern=".faa")
taxa <- taxa[-which(taxa == "Ilysp1_GeneCatalog_proteins_20121116.faa")]

protein.lengths <- list()

#ADD READING IN NUCLEOTIDES AS WELL AS PROTEINS
for (i in taxa) {
  
  print(i)
  
  proteins <- read.fasta(paste0("orthology_inference/", i), seqtype="AA")
  names(proteins) <- sub(".*\\.faa_", "", names(proteins))
  
  #Create bar to show progress
  progress.bar <- txtProgressBar(1, length(proteins), initial=0, char="=", style=3)
  
  counter <- 0
  
  for (j in names(proteins)) {
    
    counter <- counter + 1
    
    #Update progress bar
    setTxtProgressBar(progress.bar, counter)
    
    protein.lengths[[i]][j] <- length(proteins[[j]])
    
  }
  
  effectors <- scan(paste0("effector_prediction/", i, "_candidate_effectors"), character(), quote="")
  assign(paste0(i, "effectors"), effectors)
  
}

lengths.df <- stack(protein.lengths)

rm(effectors)

effectors <- unlist(mget(ls(pattern="effectors")))

lengths.df$group[!is.na(match(rownames(lengths.df), effectors))] <- "effector"
lengths.df$group[is.na(lengths.df$group)] <- "other"

lengths.labels <- data.frame(table(lengths.df$group))
lengths.labels$mean <- t.test(lengths.df$values[lengths.df$group == "effector"], lengths.df$values[lengths.df$group == "other"])$estimate
lengths.labels$label <- paste0("n=", prettyNum(lengths.labels$Freq, big.mark=","))

ggplot(lengths.df, aes(x=group, y=values)) +
  geom_hline(yintercept=300, lty="dashed") +
  geom_violin() +
  geom_boxplot(width=0.2) +
  geom_text(data=lengths.labels, aes(x=Var1, y=mean, label=label), nudge_x=0.1, inherit.aes=FALSE) +
  coord_cartesian(clip="off", ylim=c(0,2000))

#RUN GC TESTS

orthogroups <- sub("_aln_nuc.fa", "", list.files("selection/alignments/codon/", pattern="_aln_nuc.fa"))

GC.content <- list()

#Create bar to show progress
progress.bar <- txtProgressBar(1, length(orthogroups), initial=0, char="=", style=3)

for (i in orthogroups) {
  
  nucleotides <- tryCatch({read.fasta(paste0("selection/alignments/codon/", i, "_aln_nuc.fa"))}, error=function(e){NULL})
  
  #Update progress bar
  setTxtProgressBar(progress.bar, which(orthogroups == i))
  if (!is.null(nucleotides)) {
    if (sum(lengths(nucleotides)) != 0) {
      for (j in names(nucleotides)) {
        
        GC.content[[j]][i] <- GC(nucleotides[[j]])
        
      }
    }
  }
}

GC.content <- GC.content[-which(names(GC.content) == "Ilysp1_GeneCatalog_proteins_20121116")]

gc.df <- stack(GC.content)

core.effectors <- as.vector(read.csv("selection/orthogroups_absrel_meme.csv", header=FALSE)$V1)

gc.df <- as.data.frame(GC.content)

gc.df$ortho <- rownames(gc.df)

gc.df <- melt(gc.df)

gc.df$group[!is.na(match(gc.df$ortho, core.effectors))] <- "effector"

t.test(gc.df$value[gc.df$group == "effector"], gc.df$value[is.na(gc.df$group)])

ggplot(gc.df, aes(x=group, y=value)) +
  #geom_hline(yintercept=300, lty="dashed") +
  geom_violin() +
  geom_boxplot()

#PARSE BUSTED RESULTS

busted.p <- list()

#Create bar to show progress
progress.bar <- txtProgressBar(1, length(orthogroups), initial=0, char="=", style=3)

for (i in orthogroups) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, which(orthogroups == i))
  
  busted.results <- tryCatch({fromJSON(paste0("selection/hyphy/busted/", i, "_BUSTED.json"))}, error=function(e){NULL})
  
  if (!is.null(busted.results)) {
    busted.p[[i]] <-busted.results$`test results`$`p-value`
  }
}



busted.df <- stack(busted.p)
busted.df$group[!is.na(match(busted.df$ind, core.effectors))] <- "effector"
busted.df$selection[which(busted.df$values <= 0.05)] <- "Y"

ggplot(busted.df, aes(x=group, fill=selection)) +
  geom_bar(stat="count")
