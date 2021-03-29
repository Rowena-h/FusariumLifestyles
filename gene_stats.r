library(seqinr)

list <- "GCA_000149555.1_ASM14955v1_protein.faa"

lengths <- list()

#ADD READING IN NUCLEOTIDES AS WELL AS PROTEINS
for (i in list) {
  
  seq <- read.fasta(i, seqtype="AA")
  names(seq) <- sub(".*\\.faa_", "", names(seq))
  
  #Create bar to show progress
  progress.bar <- txtProgressBar(1, length(seq), initial=0, char="=", style=3)
  
  counter <- 0
  
  for (j in names(seq)) {
    
    counter <- counter + 1
    
    #Update progress bar
    setTxtProgressBar(progress.bar, counter)
    
    lengths[[i]][j] <- length(seq[[j]])
    
  }
}

#ADD LOOP READING IN EFFECTOR LISTS
effectors <- scan("../Secretome/GCA_000149555.1_ASM14955v1_protein.faa_candidate_effectors", character(), quote="")

df <- as.data.frame(lengths[1])

df$group[!is.na(match(rownames(df), effectors))] <- "effector"

t.test(df$GCA_000149555.1_ASM14955v1_protein.faa[df$group == "effector"], df$GCA_000149555.1_ASM14955v1_protein.faa[is.na(df$group)])

#RUN GC TESTS

#PARSE BUSTED RESULTS


library(ggplot2)

ggplot(df, aes(x=group, y=GCA_000149555.1_ASM14955v1_protein.faa)) +
  geom_violin() +
  geom_boxplot()
