#################################################################
#################################################################
####                                                         ####
####  Script to perform statistical tests and plot figures   ####
####                                                         ####
#################################################################
#################################################################

#Written in R v4.0.4
library(ape)          #5.6-1
library(ART)          #1.0
library(cowplot)      #1.1.1
library(deeptime)     #0.0.6.0
library(dendextend)   #1.15.2
library(plyr)         #1.8.6  must be loaded before dplyr to avoid clashes
library(dplyr)        #1.0.6
library(eulerr)       #6.1.0
library(ggplot2)      #3.3.3
library(ggalluvial)   #0.12.3
library(ggforce)      #0.3.2.9000
library(ggnewscale)   #0.4.6
library(ggplotify)    #0.0.7
library(ggpubr)       #0.4.0
library(ggrepel)      #0.9.1
library(ggthemes)     #4.2.4
library(ggtree)       #2.4.2
library(grid)         #4.0.4
library(jsonlite)     #1.7.2
library(matrixStats)  #0.61.0
library(MCMCtreeR)    #1.1
library(metR)         #0.9.2
library(multcompView) #0.1-8
library(nlme)         #3.1-152
library(pBrackets)    #1.0.1
library(phangorn)     #2.7.0
library(phytools)     #0.7-80
library(reshape2)     #1.4.4
library(rstatix)      #0.7.0
library(scales)       #1.1.1
library(seqinr)       #4.2-8
library(stringi)      #1.6.2
library(stringr)      #1.4.0
library(tidyr)        #1.1.3
library(vegan)        #2.5-7

#Colour palette
show_col(colorblind_pal()(8))

#Read in orthogroup data
load("CSEP_CAZyme_prediction/orthogroup-matrices-2022-02-10.RData")

#Read in sample metadata
metadata <- read.csv("metadata.csv")

#Make dataframe of lifestyle colours
col.df <- data.frame(lifestyle=c("endophyte", "animal pathogen", "human pathogen","animal associate",
                                 "insect mutualist", "plant associate", "plant pathogen", "saprotroph",
                                 "mycoparasite"),
                     colour=c("#009E73", "#FFE983", "#000000", "#F1BCF4", "#56B4E9",
                              "#9AE324", "dimgrey", "#0072B2", "#D55E00"))


############################    FIGURE 1    ####################################
#####################   TOPOLOGY COMPARISON   ##############################

#Read in trees
astral <- read.tree("phylogenomics/species_tree/astral/fus_proteins_62T_astral.tre")
astral$edge.length <- rep(1, length(astral$edge.length))
raxmlng <- read.tree("phylogenomics/species_tree/raxml-ng/fus_proteins_62T.raxml.support")
iqtree <- read.tree("phylogenomics/species_tree/iqtree/fus_proteins_62T_iqtree.contree")
astral.pro <- read.tree("phylogenomics/species_tree/astral/fus_proteins_62T_astralpro_multicopy.tre")
astral.pro$tip.label <- gsub("-", "_", astral.pro$tip.label)
astral.pro$tip.label[which(astral.pro$tip.label == "Ilysp1_GeneCatalog_proteins")] <- "Ilysp1_GeneCatalog_proteins_20121116"
astral.pro$tip.label <- metadata$file2[match(astral.pro$tip.label, metadata$tip)]
astral.pro$edge.length <- rep(1, length(astral.pro$edge.length))
stag <- read.tree("orthology_inference/OrthoFinder/Results_Oct22/Species_Tree/SpeciesTree_rooted.txt")

#Make vector with outgroup
outgroup <- "Ilyonectria sp."

#Create groupings for whether taxa are Fusarium sensu stricto or allied genera
allied.group <- list(allied=metadata$name[-union(setdiff(grep("Fusarium", metadata$name),
                                                         grep("'", metadata$name)), which(metadata$ingroup == 0))],
                     fusarium=metadata$name[union(setdiff(grep("Fusarium", metadata$name),
                                                          grep("'", metadata$name)), which(metadata$ingroup == 0))])

#For each species-tree method...
for (i in c("iqtree", "raxmlng", "astral", "astral.pro", "stag")) {
  
  #Get the tree
  tree <- get(i)
  tree$tip.label <- metadata$name[match(tree$tip.label, metadata$file2)]
  #Root tree
  tree <- root(tree, outgroup, resolve.root=TRUE, edgelabel=TRUE)
  #Ultrametrise
  chrono <- suppressWarnings(chronos(tree))
  #Convert to dendrogram for tanglegrams
  dend <- as.dendrogram(chrono)
  labels(dend) <- metadata$tiplab2[match(labels(dend), metadata$name)]
  
  assign(paste0(i, ".dend"), dend)
  
  #Plot tree
  gg <- ggtree(groupOTU(tree, allied.group), aes(colour=group), branch.length="none") %<+% metadata
  #Capture data structure
  gg.tree.data <- gg[["data"]] %>%
    arrange(y)
  
  #Make dataframe of species complex nodes
  sc.df <- data.frame(sc=unique(metadata$speciescomplex.abb),
                      node=NA)
  
  #Get nodes for each species complex
  for (j in 1:length(sc.df$sc)) {
    sc.df$node[j] <- MRCA(tree, metadata$name[metadata$speciescomplex.abb == sc.df$sc[j]])
  }
  
  #Collapse and scale nodes for summary tree
  gg.summary <- gg
  
  for (k in 1:length(sc.df$sc)) {
    
    if (table(gg.tree.data$speciescomplex.abb)[
      which(names(table(gg.tree.data$speciescomplex.abb)) == sc.df$sc[k])] > 1) {
      
      gg.summary <- scaleClade(gg.summary,
                               node=sc.df$node[k],
                               scale=0.4)
    }
    
  }
  
  for (k in 1:length(sc.df$sc)) {
    
    if (table(gg.tree.data$speciescomplex.abb)[
      which(names(table(gg.tree.data$speciescomplex.abb)) == sc.df$sc[k])] > 1) {
      
      gg.summary <- ggtree::collapse(gg.summary,
                                     node=sc.df$node[k],
                                     mode="max",
                                     clade_name=sc.df$sc[k],
                                     colour="black",
                                     fill="white")
      
    }
    
  }
  
  #Plot summary tree
  gg.summary <- gg.summary +
    xlim(0, 25) +
    geom_cladelab(data=sc.df,
                  mapping=aes(node=node,
                              label=sc),
                  barsize=0,
                  fontsize=1.5,
                  align=TRUE,
                  fontface="bold") +
    coord_cartesian(clip="off") +
    scale_colour_manual(values=c("red", "black")) +
    theme(legend.position="none")
  
  assign(paste0(i, ".tree"), tree)
  assign(paste0("gg.", i), gg)
  assign(paste0("gg.summary.", i), gg.summary)
  assign(paste0("gg.tree.data.", i), gg.tree.data)
  
}

#Calculate topological distances (RF) between all trees
tree.comp.mat <- signif(as.matrix(RF.dist(c(raxmlng.tree, iqtree.tree, astral.tree, astral.pro.tree, stag.tree),
                                          normalize=TRUE)), 1)

rownames(tree.comp.mat) <- c("RAxML-NG", "IQTREE", "ASTRAL-III", "ASTRAL-Pro", "STAG")
colnames(tree.comp.mat) <- c("RAxML-NG", "IQTREE", "ASTRAL-III", "ASTRAL-Pro", "STAG")

#Convert to dataframe for plotting
tree.comp.df <- melt(tree.comp.mat)
tree.comp.df <- tree.comp.df[!duplicated(data.frame(t(apply(tree.comp.df, 1, sort)))),]
tree.comp.df <- tree.comp.df[tree.comp.df$value != 0,]

#At worst, number of bipartitions which trees differed in
print(paste0(max(tree.comp.df$value) * length(bitsplits(unroot(raxmlng.tree))[3]$freq), "/",
       length(bitsplits(unroot(raxmlng.tree))[3]$freq)))

#Make dataframe to label plot
tree.comp.labels <- data.frame(method=c("Concatenated", "Coalescent", "Single-copy", "Multi-copy"), 
                               vert.x=c(0.3, 0.3, 0.1, 0.1),
                               vert.xend=c(0.3, 0.3, 0.1, 0.1),
                               vert.y=c(0.5, 1.5, 0.5, 2.5),
                               vert.yend=c(1.5, 4.5, 2.5, 4.5),
                               hor.x=c(0.5, 2.5, 0.5, 3.5),
                               hor.xend=c(2.5, 4.5, 3.5, 4.5),
                               hor.y=c(4.7, 4.7, 4.9, 4.9),
                               hor.yend=c(4.7, 4.7, 4.9, 4.9))

#Plot grid
gg.tree.comp <- ggplot(tree.comp.df, aes(Var2, Var1)) +
  geom_tile(aes(fill=value), colour="dimgrey", size=1, alpha=0.5, show.legend=FALSE) +
  geom_text(aes(label=value), size=1.5, show.legend=FALSE) +
  geom_segment(data=tree.comp.labels,
               aes(x=vert.x, xend=vert.xend, y=vert.y, yend=vert.yend, colour=method),
               size=1.5,
               inherit.aes=FALSE) +
  geom_segment(data=tree.comp.labels,
               aes(x=hor.x, xend=hor.xend, y=hor.y, yend=hor.yend, colour=method),
               size=1.5,
               inherit.aes=FALSE) +
  scale_x_discrete(position="top") +
  scale_fill_gradient(low="white", high="dimgrey") +
  scale_colour_manual(values=c("#E69F00", "#009E73", "#56B4E9", "#F0E442")) +
  theme_minimal() + 
  theme(legend.position=c(0.8, 0.25),
        legend.title=element_blank(),
        legend.key.height=unit(0.2, "cm"),
        legend.text=element_text(size=4.5, margin=margin(0, 0, 0, -5)),
        axis.text.x.top=element_text(size=3.5, face="bold", margin=margin(0, 0, 5, 0)),
        axis.text.y=element_text(face="bold", angle=90, hjust=0.5, vjust=0,
                                 size=3.5, margin=margin(0, 5, 0, 0)),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid=element_blank()) +
  coord_fixed(clip="off")

#Write topology comparison to file (Figure 1)
#tiff(file=paste0("Fig1-", Sys.Date(), ".tiff"),
#     height=2, width=6.75, unit="in", res=600, compression="lzw")

plot_grid(gg.tree.comp,
          plot_grid(gg.summary.raxmlng +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=7,
                                       yend=7),
                                   colour="#009E73",
                                   size=1.5) +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=6,
                                       yend=6),
                                   colour="#F0E442",
                                   size=1.5) +
                      annotate("text", label="RAxML-NG", size=2, angle=90, fontface="bold", x=0.5, y=32),
                    gg.summary.iqtree +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=7,
                                       yend=7),
                                   colour="#009E73",
                                   size=1.5) +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=6,
                                       yend=6),
                                   colour="#F0E442",
                                   size=1.5) +
                      annotate("text", label="IQ-TREE", size=2, angle=90, fontface="bold", x=0.5, y=32),
                    gg.summary.astral +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=7,
                                       yend=7),
                                   colour="#E69F00",
                                   size=1.5) +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=6,
                                       yend=6),
                                   colour="#F0E442",
                                   size=1.5) +
                      annotate("text", label="ASTRAL-III", size=2, angle=90, fontface="bold", x=0.5, y=32),
                    gg.summary.astral.pro +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=7,
                                       yend=7),
                                   colour="#E69F00",
                                   size=1.5) +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=6,
                                       yend=6),
                                   colour="#56B4E9",
                                   size=1.5) +
                      annotate("text", label="ASTRAL-Pro", size=2, angle=90, fontface="bold", x=0.5, y=32),
                    gg.summary.stag +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=7,
                                       yend=7),
                                   colour="#E69F00",
                                   size=1.5) +
                      geom_segment(aes(x=-Inf, xend=Inf,
                                       y=6,
                                       yend=6),
                                   colour="#56B4E9",
                                   size=1.5) +
                      annotate("text", label="STAG", size=2, angle=90, fontface="bold", x=0.5, y=32),
                    nrow=1),
          nrow=1, labels="AUTO", label_size=10, rel_widths=c(1, 2.8))

#dev.off()


## Comparison of trimming tools

raxmlng.bmge <- read.tree("phylogenomics/species_tree/raxml-ng/fus_proteins_bmge_62T.raxml.support")
iqtree.bmge <- read.tree("phylogenomics/species_tree/iqtree/fus_proteins_bmge_62T_iqtree.contree")
astral.bmge <- read.tree("phylogenomics/species_tree/astral/fus_proteins_bmge_62T_astral.tre")
astral.bmge$edge.length <- rep(1, length(astral.bmge$edge.length))

for (i in c("iqtree", "raxmlng", "astral")) {
  
  #Get the tree
  tree <- get(paste0(i, ".bmge"))
  #Edit tip names
  tree$tip.label <- metadata$name[match(tree$tip.label, metadata$tip)]
  #Root tree
  tree <- root(tree, outgroup, resolve.root=TRUE, edgelabel=TRUE)
  #Ultrametrise
  chrono <- suppressWarnings(chronos(tree))
  #Convert to dendrogram for tanglegram
  dend <- as.dendrogram(chrono)
  labels(dend) <- metadata$tiplab2[match(labels(dend), metadata$name)]
  #Untangle dendrograms from both trimming tools
  dend.comp <- untangle_labels(get(paste0(i, ".dend")), dend, method="step2side")
  
  #Plot first tree
  tr1 <- ggtree(as.phylo(dend.comp[[1]]), branch.length = "none",
                aes(lty=node %in% 
                      matchNodes(
                        as.phylo(dend.comp[[1]]),
                        as.phylo(dend.comp[[2]]))[
                          which(is.na(matchNodes(
                            as.phylo(dend.comp[[1]]),
                            as.phylo(dend.comp[[2]]))[,2])),1])) +
    xlim(0, 35) +
    geom_tiplab(size=2) +
    scale_y_continuous(expand=c(0, 2)) +
    theme_void() +
    theme(legend.position="none")
  
  #Plot second tree
  tr2 <- ggtree(as.phylo(dend.comp[[2]]), branch.length = "none",
                aes(lty=node %in% 
                      matchNodes(
                        as.phylo(dend.comp[[2]]),
                        as.phylo(dend.comp[[1]]))[
                          which(is.na(matchNodes(
                            as.phylo(dend.comp[[2]]),
                            as.phylo(dend.comp[[1]]))[,2])),1])) +
    geom_tiplab(size=2, hjust=1) +
    scale_x_reverse() +
    xlim(35, 0) +
    scale_y_continuous(expand=c(0, 2)) +
    theme_void() +
    theme(legend.position="none")
  
  #Plot dataframe for tanglegram lines
  lines.df <- data.frame(y1=tr1$data$y[tr1$data$isTip == TRUE],
                         y2=tr2$data$y[match(tr1$data$label[tr1$data$isTip == TRUE],
                                             tr2$data$label[tr2$data$isTip == TRUE])],
                         col="same")
  lines.df$col[which(lines.df$y1 != lines.df$y2)] <- "diff"
  
  if ("diff" %in% lines.df$col) {
    
    #Plot lines
    gg.lines <- ggplot() +
      geom_segment(data=lines.df, 
                   aes(x=0,
                       y=y1,
                       xend=1,
                       yend=y2,
                       col=col)) +
      scale_colour_manual(values=c("black", "grey")) +
      scale_y_continuous(expand=c(0, 2)) +
      theme_void() +
      theme(legend.position="none")
    
  } else {
    
    gg.lines <- ggplot() +
      geom_segment(data=lines.df, 
                   aes(x=0,
                       y=y1,
                       xend=1,
                       yend=y2),
                   col="grey") +
      scale_y_continuous(expand=c(0, 2)) +
      theme_void() +
      theme(legend.position="none")
    
  }
  
  #Combine into tanglegram
  tanglegram <- plot_grid(tr1, gg.lines, tr2, ncol=3, rel_widths=c(2, 0.3, 2), align="h", axis="bt")
  
  assign(paste0(i, ".trim.tanglegram"), tanglegram)
  
}

#Write tanglegram to file (Supplementary Figure 2)
#tiff(file=paste0("SupplementaryFig2-", Sys.Date(), ".tiff"),
#     height=15, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(labels="AUTO", label_size=10, ncol=1, rel_heights=c(1,1,1.1),
          raxmlng.trim.tanglegram, iqtree.trim.tanglegram,
          ggdraw(add_sub(astral.trim.tanglegram, "TrimAl\t\t\t\t\t\t\t\t\t\t\t\t\tBMGE", size=10)))
#dev.off()

####################################

#Check distributions of branch lengths
branch.lengths.df <- data.frame(raxmlng=raxmlng.tree$edge.length, iqtree=iqtree.tree$edge.length, stag=stag.tree$edge.length)

branch.lengths.df2 <- melt(branch.lengths.df)

ggplot(branch.lengths.df2, aes(x=value)) +
  facet_wrap(~ variable,
             labeller=labeller(variable=c(raxmlng="RAxML-NG", iqtree="IQ-TREE", stag="STAG"))) +
  geom_density(alpha=0.5, fill="grey") +
  labs(x="Branch length (#substitutions/site)", y="Density") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text=element_blank())

###############################

#Proportion of supported branches
#RAxML-NG Felsenstein's bootstrap >= 70
suppressWarnings(round(length(which(as.numeric(raxmlng$node.label) >= 70)) /
                         length(na.omit(as.numeric(raxmlng$node.label))) * 100))
#IQ-TREE ultrafast bootstrap >= 95
suppressWarnings(round(length(which(as.numeric(iqtree$node.label) >= 95)) /
                         length(na.omit(as.numeric(iqtree$node.label))) * 100))
#ASTRAL-III LPP >= 0.95
suppressWarnings(round(length(which(as.numeric(astral$node.label) >= 0.95)) /
                         length(na.omit(as.numeric(astral$node.label))) * 100))
#ASTRAL-Pro LPP >=0.95
suppressWarnings(round(length(which(as.numeric(astral.pro$node.label) >= 0.95)) /
                         length(na.omit(as.numeric(astral.pro$node.label))) * 100))
#STAG proportion of trees with same branch >=0.30
suppressWarnings(round(length(which(as.numeric(stag$node.label) >= 0.30)) /
                         length(na.omit(as.numeric(stag$node.label))) * 100))


############################    FIGURE 2    ####################################
#####################   PHYLOGENY AND LIFESTYLE   ##############################

## MCMCTREE CONVERGENCE PLOTS ##

#For each relax clock model...
for (i in c("independent", "correlated")) {
  
  #Read in the generation data
  mcmc1.gens <- read.csv(paste0("divergence_time_estimation/mcmctree/run1_", i, "/mcmc_run1_", i, ".txt"), sep="\t")
  mcmc2.gens <- read.csv(paste0("divergence_time_estimation/mcmctree/run2_", i, "/mcmc_run2_", i, ".txt"), sep="\t")
  
  #Calculate the effective sample size
  ESS <- mean(apply(mcmc1.gens[,-1], 2, effectiveSize))
  
  #Make dataframe of posterior means for both chains
  mcmc.df <- data.frame(run1=colMeans(mcmc1.gens[2:62]),
                        run2=colMeans(mcmc2.gens[2:62]))
  
  #Plot chain convergence
  gg.mcmc <- ggplot(mcmc.df, aes(x=run1, y=run2)) +
    geom_abline(colour="dimgrey") +
    geom_point() +
    labs(x="Posterior mean chain 1 (100MY)",
         y="Posterior mean chain 2 (100MY)",
         subtitle=paste0("ESS=", round(ESS))) +
    coord_fixed() +
    theme(axis.title=element_text(size=6),
          axis.text=element_text(size=5),
          plot.subtitle=element_text(size=8))
  
  assign(paste0("gg.mcmc.", i), gg.mcmc)
  
  #Read in and format MCMCTRee output
  mcmc1.out <- read.csv(paste0("divergence_time_estimation/mcmctree/run1_", i, "/mcmctree_step2_out.txt"),
                        sep="\t", skip=616, nrows=61, header=FALSE)
  mcmc1.out <- as.data.frame(do.call(rbind, strsplit(mcmc1.out$V1, "\\s+")))
  mcmc1.out <- mcmc1.out[1:7]
  colnames(mcmc1.out) <- c("node", "posterior_mean", "95_equal_tail_CI_lower", "95_equal_tail_CI_upper",
                           "95_HPD_CI_lower", "95_HPD_CI_upper", "HPD_CI_width")
  mcmc1.out[] <- lapply(mcmc1.out, function(x) gsub("\\(", "", x))
  mcmc1.out[] <- lapply(mcmc1.out, function(x) gsub("\\)", "", x))
  mcmc1.out[] <- lapply(mcmc1.out, function(x) gsub(",", "", x))
  mcmc1.out[2:7] <- data.frame(apply(mcmc1.out[2:7], 2, function(x) as.numeric(as.character(x))))
  
  #Calculate confidence interval width
  mcmc1.out$HPD_CI_width <- mcmc1.out$`95_HPD_CI_upper` - mcmc1.out$`95_HPD_CI_lower`
  
  #Plot infinite-sites plot
  gg.infinitesite <- ggplot(mcmc1.out, aes(x=posterior_mean, y=HPD_CI_width)) +
    geom_smooth(method="lm", se=FALSE, colour="dimgrey") +
    geom_point() +
    labs(x="Posterior mean chain 1 (100MY)",
         y="Confidence interval width (100MY)") +
    theme(axis.title=element_text(size=6),
          axis.text=element_text(size=5))
  
  assign(paste0("gg.infinitesite.", i), gg.infinitesite)
  
}

#Write to file (Supplementary Figure 12)
#tiff(file=paste0("SupplementaryFig12-", Sys.Date(), ".tiff"),
#     height=4, width=6.75, units="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.mcmc.correlated, gg.infinitesite.correlated, align="h", axis="tb"),
          plot_grid(gg.mcmc.independent, gg.infinitesite.independent, align="h", axis="tb"), nrow=2,
          labels="AUTO", label_size=10)
#dev.off()


## DATED PHYLOGENY ##

#Plot species tree
gg.speciestree <- ggtree(raxmlng.tree, branch.length="none") %<+% metadata

#Capture species tree data structure
gg.tree.data.speciestree <- gg.speciestree[["data"]] %>%
  arrange(y)

#For analyses from both relaxed molecular clock methods...
for (clock in c("independent", "correlated")) {
  
  #Read in MCMCTree tree
  dated.tree <- readMCMCtree(paste0("divergence_time_estimation/mcmctree/run1_", clock, "/FigTree.tre"),
                             forceUltrametric=TRUE)
  
  #Convert branch length unit to 1MY
  dated.tree$apePhy$edge.length <- dated.tree$apePhy$edge.length * 100
  #Replace tip labels with names
  dated.tree$apePhy$tip.label <- metadata$name[match(dated.tree$apePhy$tip.label, metadata$short.tip)]
  #Add ML support values from to correct nodes
  nodes <- matchNodes(raxmlng.tree, dated.tree$apePhy)
  dated.tree$apePhy$node.label <- gg.tree.data.speciestree$label[
    match(nodes[,1], gg.tree.data.speciestree$node)][order(nodes[,2])]
  
  assign(paste0("dated.tree.", clock), dated.tree)
  
  #Plot tree
  gg.datedtree <- ggtree(dated.tree$apePhy, linetype=NA) %<+% metadata
  
  #Create vector with the order of the tree tips for plotting
  gg.tree.data.dated <- gg.datedtree[["data"]] %>%
    arrange(y)
  tip.order.datedtree <- rev(gg.tree.data.dated$label[gg.tree.data.dated$isTip == "TRUE"])
  
  #Create vector of fontface for new genomes in this study
  label.face.datedtree <- rev(metadata$new[match(tip.order.datedtree, metadata$name)])
  
  for (i in 1:length(label.face.datedtree)){
    if (label.face.datedtree[i] == "Y") {
      label.face.datedtree[i] <- "bold.italic"
    } else {
      label.face.datedtree[i] <- "italic"
    }
  }
  
  tiplabel.face.datedtree <- rep("italic", length(dated.tree$apePhy$tip.label))
  tiplabel.face.datedtree[which(metadata$new[
    match(dated.tree$apePhy$tip.label, metadata$name)] == "Y")] <- "bold.italic"
  
  #Make dataframe for plotting divergence time confidence intervals
  confidence.intervals <- as.data.frame(dated.tree$nodeAges)
  confidence.intervals <- confidence.intervals * 100
  confidence.intervals$ymin <- gg.tree.data.dated$y[match(rownames(confidence.intervals),
                                                          gg.tree.data.dated$node)] -0.2
  confidence.intervals$ymax <- gg.tree.data.dated$y[match(rownames(confidence.intervals),
                                                          gg.tree.data.dated$node)] +0.2
  
  #Add time scale to tree plot
  gg.datedtree <- gg.datedtree +
    coord_geo(clip="off",
              xlim=c(max(confidence.intervals$`95%_upper`) * -1, 70),
              ylim=c(0, Ntip(dated.tree$apePhy)),
              dat=list("epochs", "periods"),
              pos=list("bottom", "bottom"),
              size=list(1,1.5),
              height=list(unit(0.25, "line"), unit(0.5, "line")),
              abbrv=list(TRUE, FALSE),
              center_end_labels=TRUE,
              lwd=0, alpha=0.5, neg=TRUE, expand=FALSE) +
    scale_x_continuous(breaks=rev(seq(0, round(max(confidence.intervals$`95%_upper`)), 10)) * -1,
                       labels=rev(seq(0, round(max(confidence.intervals$`95%_upper`)), 10)),
                       expand=c(0, 0),
                       name="Million years") +
    theme(axis.text.x.bottom=element_text(size=5),
          axis.title.x.bottom=element_text(size=5))
  
  #Reverse the x axis
  gg.datedtree <- revts(gg.datedtree)
  
  #Make dataframe of species complex nodes
  sc.df.dated <- data.frame(sc=unique(metadata$speciescomplex.abb),
                            node=NA)
  
  #Get nodes for each species complex
  for (i in 1:length(sc.df.dated$sc)) {
    sc.df.dated$node[i] <- MRCA(dated.tree$apePhy,
                                   metadata$name[metadata$speciescomplex.abb == sc.df.dated$sc[i]])
  }
  
  #Make alternated coding for highlights on tree
  sc.df.dated <- sc.df.dated[match(na.omit(unique(gg.tree.data.dated$speciescomplex.abb)), sc.df.dated$sc),]
  sc.df.dated$box <- rep(c(0,1), length.out=length(sc.df.dated$sc))
  
  #Add final annotations to tree plot
  gg.datedtree <- gg.datedtree +
    geom_highlight(data=sc.df.dated, 
                   aes(node=node, fill=as.factor(box)), alpha=1, extend=200, show.legend=FALSE) +
    geom_cladelab(data=sc.df.dated,
                  mapping=aes(node=node, label=sc),
                  offset.text=2, offset=70, align=TRUE, barsize=0.2, fontsize=1.3) +
    geom_rect(data=confidence.intervals,
              aes(xmin=`95%_lower` * -1, 
                  xmax=`95%_upper` * -1,
                  ymin=ymin,
                  ymax=ymax),
              alpha=0.1,
              size=0.1,
              colour="grey",
              inherit.aes=FALSE) +
    geom_point(data=confidence.intervals,
               aes(x=mean * -1,
                   y=ymin + 0.2),
               colour="grey",
               size=0.5,
               inherit.aes=FALSE) +
    scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
    geom_vline(xintercept=epochs$max_age[which((epochs$max_age) < max(confidence.intervals$`95%_upper`))] * -1,
               linetype="dashed", colour="grey", size=0.2) +
    geom_tiplab(aes(label=tiplab), fontface=tiplabel.face.datedtree, size=1.5, offset=1) +
    geom_tippoint(aes(colour=lifestyle), size=0.7, show.legend=FALSE) +
    scale_colour_manual(values=col.df$colour[match(sort(unique(metadata$lifestyle[match(dated.tree$apePhy$tip.label,
                                                                                        metadata$name)])),
                                                   col.df$lifestyle)]) +
    new_scale_colour() +
    geom_tree(size=0.25, aes(colour=node %in% gg.datedtree$data$node[
      suppressWarnings(intersect(which(as.numeric(gg.datedtree$data$label) < 70),
                                 which(gg.datedtree$data$isTip == FALSE)))]),
      show.legend=FALSE) +
    scale_colour_manual(values=c("black", "red"))
  
  #Format confidence intervals for table
  confidence.intervals.tab <- data.frame(Node=as.numeric(rownames(confidence.intervals)) - Ntip(dated.tree$apePhy),
                                         Mean.age=round(confidence.intervals$mean, digits=1),
                                         `95.HPD`=paste0(round(confidence.intervals$`95%_lower`, digits=1), ", ",
                                                         round(confidence.intervals$`95%_upper`, digits=1)))
  
  confidence.intervals.tab <- split(confidence.intervals.tab,
                                    rep(c("first", "second"),
                                        each=ceiling(length(confidence.intervals.tab$Node)/2)))
  
  #Plot confidence intervals table in two separate columns
  for (i in c("first", "second")) {
    
    gg.ci <- ggtexttable(confidence.intervals.tab[[i]],
                         rows=NULL,
                         theme=ttheme(base_size=3,
                                      padding=unit(c(1, 1), "mm")),
                         cols=c("Node", "Mean age", "95% HPD"))
    
    assign(paste0("gg.ci.", i, ".", clock), gg.ci)
    
  }
  
  assign(paste0("gg.datedtree.", clock), gg.datedtree)
  assign(paste0("gg.tree.data.dated.", clock), gg.tree.data.dated)
  
}

#Write dated trees with confidence interval tables to file (Supplementary Figure 3)
#tiff(file=paste0("SupplementaryFig3-", Sys.Date(), ".tiff"),
#     height=8, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.datedtree.correlated +
                      geom_text_repel(data=gg.tree.data.dated.correlated[
                        gg.tree.data.dated.correlated$node %in% c(64, 66),],
                        aes(x=x-max(gg.tree.data.dated.correlated$x), y=y,
                            label=c("Generic concept\nsensu lato",
                                    "Generic concept\nsensu stricto")),
                        hjust=1,
                        size=1.1,
                        segment.size=0.25,
                        min.segment.length=unit(0, 'lines'),
                        xlim=c(-Inf, Inf),
                        nudge_x=-3,
                        nudge_y=2.5) +
                      geom_label2(aes(subset=!isTip, label=node - Ntip(dated.tree$apePhy)),
                                  size=1, label.padding=unit(1, "pt")) +
                      theme(plot.margin=margin(0.8, 1, 0, 0, unit="cm")),
                    gg.ci.first.correlated, gg.ci.second.correlated,
                    rel_widths=c(6, 1 ,1), nrow=1, align="h", axis="tb"),
          plot_grid(gg.datedtree.independent +
                      geom_text_repel(data=gg.tree.data.dated.independent[
                        gg.tree.data.dated.independent$node %in% c(64, 66),],
                        aes(x=x-max(gg.tree.data.dated.independent$x), y=y,
                            label=c("Generic concept\nsensu lato",
                                    "Generic concept\nsensu stricto")),
                        hjust=1,
                        size=1.1,
                        segment.size=0.25,
                        min.segment.length=unit(0, 'lines'),
                        xlim=c(-Inf, Inf),
                        nudge_x=-3,
                        nudge_y=2.5) +
                      geom_label2(aes(subset=!isTip, label=node - Ntip(dated.tree$apePhy)),
                                  size=1, label.padding=unit(1, "pt")) +
                      theme(plot.margin=margin(0.8, 1, 0, 0, unit="cm")),
                    gg.ci.first.independent, gg.ci.second.independent,
                    rel_widths=c(6, 1 ,1), nrow=1, align="h", axis="tb"),
          ncol=1, labels="AUTO", label_size=10, align="v", axis="lr")
#dev.off()


## NUMBER OF GENES/CSEPs BARGRAPH ##

#Filter orthogroup presence-absence matrix for orthogroups with at least one predicted CSEP/CAZyme
CSEP.orthogroups <- orthogroups.count.ingroup1[
  match(rownames(CSEP.count.ingroup1)[which(rowSums(CSEP.count.ingroup1) > 0)],
        rownames(orthogroups.count.ingroup1)),]
CAZyme.orthogroups <- orthogroups.count.ingroup1[
  match(rownames(CAZyme.count.ingroup1)[which(rowSums(CAZyme.count.ingroup1) > 0)],
        rownames(orthogroups.count.ingroup1)),]

#Read in phylogenetic PCA from lifestyle test results
for (i in c("CSEPs", "CAZymes", "orthogroups")) {
  
  phy.pca.result <- read.csv(paste0("lifestyle_comparison/", i, "/metadata.csv"))
  rownames(phy.pca.result) <- metadata$file2[match(phy.pca.result$genome, metadata$short.tip)]
  
  assign(paste0("phy.pca.result.", i), phy.pca.result)
  
}

#For each category...
for (i in c("orthogroup", "CSEP", "CAZyme")) {
  
  if (i == "orthogroup") {
    
    tmp.orthogroups <- orthogroups.count.ingroup1
    
  } else {
    
    tmp.orthogroups <- get(paste0(i, ".orthogroups"))
    
  }
  
  #Check for correlation between N50 and number of CSEPs/CAZymes
  tmp.df <- data.frame(num=colSums(tmp.orthogroups > 0))
  tmp.df$taxon <- rownames(tmp.df)
  tmp.df$lifestyle <- metadata$lifestyle[match(tmp.df$taxon, metadata$file2)]
  tmp.df$N50 <- metadata$N50[match(tmp.df$taxon, metadata$file2)]
  cor <- cor.test(tmp.df$num, tmp.df$N50)
  
  if (cor$p.value >= 0.05) {
    
    print(paste0("No correlation between N50 and number of ", i, "s"))
    
  } else {
    
    print(paste0("Correlation between N50 and number of ", i, "s"))
    
  }
  
  #Add PC1 and PC2
  tmp.df <- cbind(tmp.df,
                  phy.pca.result.orthogroups[match(tmp.df$taxon,
                                                   rownames(phy.pca.result.orthogroups)),
                                             c("PC1", "PC2")])
  
  #Test for significant difference in number of CSEPs/CAZymes between lifestyles
  #Check for normality of residuals
  plot(ggqqplot(residuals(lm(num ~ lifestyle, data=tmp.df))) +
         ggtitle(i))
  
  #Check for homogeneity of variances
  num.genes.levene <- tmp.df %>% levene_test(num ~ lifestyle)
  num.genes.levene <- data.frame(formula=paste0("Number of ", i, "s ~ lifestyle"), num.genes.levene)
  assign(paste0("num.genes.levene.", i), num.genes.levene)
  
  #If homogeneity of variance assumption met...
  if (num.genes.levene$p >= 0.05) {
    
    print("Homogeneity of variance assumption met, doing ANOVA")
    #Do ANOVA including lifestyle and first 2 PCs
    num.genes.anova <- tmp.df %>% anova_test(num ~ PC1 + PC2 + lifestyle)
    #Format
    num.genes.anova[, c(3, 6, 7)] <- NULL
    colnames(num.genes.anova) <- c("Effect", "Df", "F", "p")
    num.genes.anova$formula <- paste0("Number of ", i, "s ~ PC1 + PC2 + lifestyle")
    assign(paste0("num.genes.anova.", i), num.genes.anova)
    num.genes.p <- num.genes.anova$p[3]
    
    #If ANOVA is significant...
    if (num.genes.p < 0.05) {
      
      print("Number of ", i, "s significant according to ANOVA, doing TukeyHSD")
      #Do TukeyHSD pairwise tests
      num.genes.pw.df <- tmp.df %>% tukey_hsd(num ~ lifestyle)
      num.genes.pw.df[c(1, 4, 9)] <- NULL
      num.genes.pw.df$formula <- paste0("Number of ", i, "s ~ lifestyle")
      assign(paste0("num.genes.pw.df.", i), num.genes.pw.df)
      
    } else {
      
      print(paste0("ANOVA, number of ", i, "s not significant: p=", signif(num.genes.p, 1)))
      
    }
    
  } else {
    
    print("Homogeneity of variance assumption not met, doing aligned rank tranform ANOVA")
    #Do ART
    num.genes.anova <- aligned.rank.transform(num ~ PC1 + PC2 + lifestyle, data=tmp.df)
    #Format
    num.genes.anova$significance$`Sum Sq` <- NULL
    colnames(num.genes.anova$significance) <- c("Df", "F", "p")
    num.genes.anova$significance$Effect <- rownames(num.genes.anova$significance)
    num.genes.anova$significance$formula <- paste0("Number of ", i, "s ~ PC1 + PC2 + lifestyle")
    assign(paste0("num.genes.anova.", i), num.genes.anova)
    num.genes.p <- num.genes.anova$significance$p[3]
    
    #If ANOVA is significant...
    if (num.genes.p < 0.05) {
      
      print("Number of ", i, "s significant according to aligned rank transform ANOVA, doing Games Howell test")
      #Do Games Howell pairwise tests
      num.genes.pw.df <- tmp.df %>% games_howell_test(num ~ lifestyle)
      num.genes.pw.df[c(1, 8)] <- NULL
      num.genes.pw.df$formula <- paste0("Number of ", i, "s ~ lifestyle")
      assign(paste0("num.genes.pw.df.", i), num.genes.pw.df)
      
    } else {
      
      print(paste0("Aligned rank transform ANOVA, number of ", i,
                   "s not significant: p=", signif(num.genes.p, 1)))
      
    }
    
  }
  
}

#Make dataframe categorising genes into core, accessory or specific
genes.df <- data.frame(taxon=rep(colnames(orthogroups.count.ingroup1), each=3), 
                       category=rep(c("specific", "accessory", "core"),
                                    length(colnames(orthogroups.count.ingroup1))),
                       orthogroups=NA,
                       CSEPs=NA,
                       CAZymes=NA)

#For each taxon...
for (i in unique(genes.df$taxon)) {
  
  #For each category...
  for (j in unique(genes.df$category)) {
    
    #Get number of orthogroups
    genes.df$orthogroups[intersect(which(genes.df$taxon == i), which(genes.df$category == j))] <-
      table(orthogroups.stats.ingroup1$category[
        match(rownames(orthogroups.count.ingroup1[orthogroups.count.ingroup1[, i] > 0,]),
              orthogroups.stats.ingroup1$orthogroup)])[j]
    
    #Get number of CSEPs
    genes.df$CSEPs[intersect(which(genes.df$taxon == i), which(genes.df$category == j))] <-
      table(orthogroups.stats.ingroup1$category[
        match(rownames(CSEP.orthogroups[CSEP.orthogroups[, i] > 0,]),
              orthogroups.stats.ingroup1$orthogroup)])[j]
    
    #Get number of CAZymes
    genes.df$CAZymes[intersect(which(genes.df$taxon == i), which(genes.df$category == j))] <-
      table(orthogroups.stats.ingroup1$category[
        match(rownames(CAZyme.orthogroups[CAZyme.orthogroups[, i] > 0,]),
              orthogroups.stats.ingroup1$orthogroup)])[j]
    
  }
}

#Add taxon name
genes.df$name <- metadata$name[match(genes.df$taxon, metadata$file2)]
#Order categories
genes.df$category <- factor(genes.df$category, levels=c("specific", "accessory", "core"))
#Order taxa by tips in the tree
genes.df$name <- factor(genes.df$name, levels=rev(tip.order.datedtree))

#Melt for plotting
genes.df <- melt(genes.df)
#Add column for order of taxa in tree
genes.df$y <- as.numeric(factor(genes.df$name))

#Make function for dynamic axis labels
addUnits <- function(n) {
  labels <- ifelse(n < 1000, n,
                   ifelse(n < 1e6, paste0(round(n/1e3), 'k')
                   )
  )
  return(labels)
}

#Plot bargraphs of number of genes
gg.gene.numbers <- ggplot(genes.df, aes(y=name, x=value, fill=category)) +
  facet_wrap(. ~ variable, scales="free",
             labeller=labeller(variable=c(CAZymes="CAZymes", CSEPs="CSEPs", orthogroups="All genes"))) +
  geom_bar(stat="identity", size=0.5, width=0.6) +
  scale_y_discrete(limits=rev(tip.order.datedtree),
                   expand=expansion(0, 0)) +
  coord_cartesian(clip="off",
                  ylim=c(0, 62)) +
  guides(fill=guide_legend(nrow=1)) +
  scale_x_continuous(expand=c(0, 0),
                     position="top",
                     labels=addUnits) +
  scale_fill_manual(values=c("lightgrey", "darkgrey", "dimgrey"),
                    breaks=c("core", "accessory", "specific"),
                    labels=c("Core", "Accessory", "Specific")) +
  theme_minimal() +
  theme(strip.placement="outside",
        strip.text=element_text(size=6, margin=margin(0, 0, 1, 0)),
        axis.title=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x.top=element_text(size=3),
        panel.grid.major.x=element_line(colour="white"),
        panel.grid.minor.x=element_line(colour="white"),
        panel.grid.major.y=element_blank(),
        legend.position=c(0.5, 0),
        legend.text=element_text(size=5, margin=margin(0, 4, 0, 0)),
        legend.key.size=unit(0.2, "cm"),
        legend.spacing.x=unit(0.1, "cm"),
        legend.title=element_blank())


## NAMED CSEPS ABUNDANCE MATRIX ##

#Proportion of CSEPs that could be functionally annotated
length(which(!is.na(orthogroups.stats.ingroup0$CSEP_name))) /
  length(which(!is.na(orthogroups.stats.ingroup0$CSEP)))

#Read in PHI-base data
phibase.df <- read.csv("CSEP_CAZyme_prediction/blastp/phi-base_current.csv")

csep.name.list <- list()

#For each taxon...
for (i in unique(genes.df$taxon)) {
  
  #Get numbers of named CSEPs
  csep.name.list[[i]] <-
    as.list(table(sub(",.*", "", orthogroups.stats.ingroup1$PHI.base_entry[
      match(rownames(orthogroups.count.ingroup1[orthogroups.count.ingroup1[, i] > 0,]),
            orthogroups.stats.ingroup1$orthogroup)])))
  
  #Format CSEP names
  names(csep.name.list[[i]]) <- sub("_.*", "", names(csep.name.list[[i]]))
  
}

#Melt into dataframe
csep.name.df <- melt(data.frame(taxon=names(csep.name.list),
                                data.table::rbindlist(csep.name.list, fill=TRUE, use.names=TRUE)))
csep.name.df$CSEP_name <- phibase.df$Gene[match(gsub("\\.", ":", csep.name.df$variable),
                                                phibase.df$PHI_MolConn_ID)]
#Add mutant phenotypes
csep.name.df$group <- phibase.df$Mutant.Phenotype[match(gsub("\\.", ":", csep.name.df$variable),
                                                        phibase.df$PHI_MolConn_ID)]

#Format group labels
csep.name.df$group <- str_to_sentence(csep.name.df$group)
csep.name.df$group <- sub("Effector \\(plant avirulence determinant\\)", "Effector", csep.name.df$group)
csep.name.df$group <- sub("Increased virulence \\(hypervirulence\\)", "Hyp", csep.name.df$group)
csep.name.df$group <- sub("Loss of pathogenicity", "Loss\npath", csep.name.df$group)
csep.name.df$group <- sub("Lethal", "L", csep.name.df$group)

#Add taxon name
csep.name.df$name <- metadata$name[match(csep.name.df$taxon, metadata$file2)]
csep.name.df$tiplab <- metadata$tiplab2[match(csep.name.df$taxon, metadata$file2)]
#Order taxa by tips in the tree
csep.name.df$name <- factor(csep.name.df$name, levels=rev(tip.order.datedtree))

#Add whether CSEP was also predicted to be a CAZyme
csep.name.df$CAZyme <- 1
for (i in unique(csep.name.df$variable)) {
  
  matches <- unique(orthogroups.stats.ingroup1$CAZyme_family[grep(gsub("\\.", ":", i),
                                                                  orthogroups.stats.ingroup1$PHI.base_entry)])
  
  if (TRUE %in% !is.na(matches)) {
    
    csep.name.df$CAZyme[csep.name.df$variable == i] <- 2
  
  }
  
}

csep.name.df$CSEP_name <- sub(" \\(.*", "", csep.name.df$CSEP_name)

#Correct or summarise gene names for plotting
csep.name.df$CSEP_name[grep("O-methylsterigmatocystin oxidoreductase", csep.name.df$CSEP_name)] <- "ordA"
csep.name.df$CSEP_name[grep("endo-1,4-beta-xylanase", csep.name.df$CSEP_name)] <- 
  gsub(".*\\[| family]", "", csep.name.df$CSEP_name[grep("endo-1,4-beta-xylanase", csep.name.df$CSEP_name)])
csep.name.df$CSEP_name <- sub("Six", "SIX", csep.name.df$CSEP_name)

#Order in plot
csep.name.df$CSEP_name <- 
  factor(csep.name.df$CSEP_name,
         levels=unique(csep.name.df$CSEP_name)[order(tolower(unique(csep.name.df$CSEP_name)))])

#Make function to bold genes which were also predicted to be CAZymes
bold_labels_cazymes <- function(breaks) {
  cazymes <- filter(csep.name.df, CSEP_name %in% breaks) %>%
    group_by(CSEP_name) %>%
    dplyr::summarise(CAZyme=mean(CAZyme)) %>%
    mutate(check=str_detect(CAZyme, "2")) %>% 
    pull(check)
  labels <- purrr::map2(
    breaks, cazymes,
    ~ if (.y) bquote(bold(.(.x))) else bquote(plain(.(.x)))
  )
  parse(text=labels)
}

#Plot as two separate grids to stack
gg.cseps.1 <- ggplot(csep.name.df[c(which(csep.name.df$group == "Effector"),
                                    which(csep.name.df$group == "Reduced virulence")),],
                     aes(x=CSEP_name, y=name, fill=value)) +
  facet_grid(. ~ group, scales="free", space="free") +
  geom_tile(color="grey", size=0.1) +
  scale_y_discrete(limits=rev(tip.order.datedtree),
                   expand=expansion(0, 0)) +
  scale_x_discrete(labels=bold_labels_cazymes) +
  coord_cartesian(clip="off",
                  ylim=c(0, 63)) +
  scale_fill_gradient(low="#F0E442", high="#CC79A7",
                      breaks=pretty_breaks(),
                      guide=guide_colourbar(title="Copynumber",
                                            title.position="left",
                                            direction="horizontal",
                                            title.hjust=0,
                                            title.vjust=0.8)) +
  theme_minimal() +
  theme(legend.position="top",
        legend.direction="horizontal",
        legend.margin=margin(5, 0, -5, 0),
        legend.title=element_text(size=5, face="bold", margin=margin(0, 10, 0, 0)),
        legend.text=element_text(size=3, margin=margin(0, 3, 0, 0)),
        legend.key.size=unit(6, "pt"),
        strip.text=element_text(size=3.5, face="bold", margin=margin(2, 0, 2, 0)),
        panel.spacing=unit(0.15, "lines"),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=3, margin=margin(-2, 0, 0, 0)))

gg.cseps.2 <- ggplot(csep.name.df[c(which(csep.name.df$group == "Unaffected pathogenicity"),
                                    which(csep.name.df$group == "Hyp"),
                                    which(csep.name.df$group == "L"),
                                    which(csep.name.df$group == "Loss\npath")),],
                     aes(x=CSEP_name, y=name, fill=value)) +
  facet_grid(. ~ group, scales="free", space="free") +
  geom_tile(color="grey", size=0.1) +
  scale_y_discrete(limits=rev(tip.order.datedtree),
                   expand=expansion(0, 0)) +
  scale_x_discrete(labels=bold_labels_cazymes) +
  coord_cartesian(clip="off",
                  ylim=c(0, 63)) +
  scale_fill_gradient(low="#F0E442", high="#CC79A7",
                      breaks=pretty_breaks()) +
  theme_minimal() +
  theme(legend.position="none",
        strip.text=element_text(size=3.5, face="bold", margin=margin(2, 0, 2, 0)),
        panel.spacing=unit(0.15, "lines"),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=3, margin=margin(-2, 0, 0, 0)))

#Plot tree for grids to be beside
gg.cseps.tree <- ggtree(dated.tree$apePhy, linetype=NA) %<+% metadata +
  geom_highlight(data=sc.df.dated, 
                 aes(node=node, fill=as.factor(box)), alpha=1, extend=500, show.legend=FALSE) +
  geom_cladelab(data=sc.df.dated,
                mapping=aes(node=node, label=sc),
                offset.text=2, offset=90, align=TRUE, barsize=0.2, fontsize=1.3) +
  coord_cartesian(clip="off") +
  scale_y_continuous(expand=expansion(0, 0),
                     limits=c(0, 63)) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tiplab(aes(label=tiplab2), fontface=tiplabel.face.datedtree, size=1.5, offset=1) +
  geom_tippoint(aes(colour=lifestyle), size=0.7, show.legend=FALSE) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(metadata$lifestyle[match(dated.tree$apePhy$tip.label,
                                                                                      metadata$name)])),
                                                 col.df$lifestyle)]) +
  new_scale_colour() +
  geom_tree(size=0.25, aes(colour=node %in% gg.datedtree$data$node[
    suppressWarnings(intersect(which(as.numeric(gg.datedtree$data$label) < 70),
                               which(gg.datedtree$data$isTip == FALSE)))]),
    show.legend=FALSE) +
  scale_colour_manual(values=c("black", "red"))

#Write to file (Supplementary Figure 6)
#tiff(file=paste0("SupplementaryFig6-", Sys.Date(), ".tiff"),
#     height=8, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.cseps.tree, gg.cseps.1,
                    nrow=1, rel_widths=c(0.8, 1), align="h", axis="bt"),
          plot_grid(gg.cseps.tree, gg.cseps.2,
                    nrow=1, rel_widths=c(0.8, 1), align="h", axis="bt"),
          nrow=2, rel_heights=c(1, 0.95)
)
#dev.off()


## CAZYMES WITH KNOWN SUBSTRATES ABUNDANCE MATRIX ##

cazyme.fam.list <- list()

#For each taxon...
for (i in unique(genes.df$taxon)) {
  
  #Get numbers of CAZymes
  cazyme.fam.list[[i]] <-
    as.list(table(sub(",.*", "", orthogroups.stats.ingroup1$CAZyme_family[
      match(rownames(orthogroups.count.ingroup1[orthogroups.count.ingroup1[, i] > 0,]),
            orthogroups.stats.ingroup1$orthogroup)])))
  
  #Format CAZyme names
  names(cazyme.fam.list[[i]]) <- gsub("\\+[[:digit:]].*\\+", "\\1\\+", names(cazyme.fam.list[[i]]))
  names(cazyme.fam.list[[i]]) <- gsub("\\+[[:digit:]].*$", "\\1\\+", names(cazyme.fam.list[[i]]))
  names(cazyme.fam.list[[i]]) <- gsub("\\+$", "", names(cazyme.fam.list[[i]]))
  
}

cazyme.fam.list <- lapply(cazyme.fam.list, function(y) {lapply(split(y, names(y)), function(x) {Reduce("+", x) })})

#Melt into dataframe
cazyme.fam.df <- melt(data.frame(taxon=names(cazyme.fam.list),
                                 data.table::rbindlist(cazyme.fam.list, fill=TRUE, use.names=TRUE),
                                 substrate=NA))
cazyme.fam.df$variable <- gsub("\\.", "+", cazyme.fam.df$variable)

#Get substrate groupings
substrates.df <- read.csv("cazyme_substrates.csv")

#Add substrates to dataframe
for (i in 1:length(unique(substrates.df$CAZy.Family))) {
  
  substrates <- substrates.df$Substrate[substrates.df$CAZy.Family == unique(substrates.df$CAZy.Family)[i]]
  
  if (length(substrates) > 0) {
    
    cazyme.fam.df$substrate[grep(paste0("\\b", unique(substrates.df$CAZy.Family)[i], "\\b"),
                                 cazyme.fam.df$variable)] <-paste(substrates, collapse=",") 
    
  }
  
}

#Split rows with more than one substrate
cazyme.fam.df <- separate_rows(cazyme.fam.df, sep=",", substrate)
#Add taxon names
cazyme.fam.df$name <- metadata$name[match(cazyme.fam.df$taxon, metadata$file2)]
cazyme.fam.df$tiplab <- metadata$tiplab2[match(cazyme.fam.df$taxon, metadata$file2)]
#Order taxa by tips in the tree
cazyme.fam.df$name <- factor(cazyme.fam.df$name, levels=rev(tip.order.datedtree))
#Abbreviate substrates for plotting
cazyme.fam.df$substrate <- sub("Cutin", "Cu", cazyme.fam.df$substrate)
cazyme.fam.df$substrate <- sub("Lignin", "L", cazyme.fam.df$substrate)

#Add whether CSEP was also predicted to be a CAZyme
cazyme.fam.df$CSEP <- 1
for (i in unique(cazyme.fam.df$variable)) {
  
  matches <- unique(orthogroups.stats.ingroup1$CSEP_name[grep(i, orthogroups.stats.ingroup1$CAZyme_family)])
  
  if (TRUE %in% !is.na(matches)) {
    
    cazyme.fam.df$CSEP[cazyme.fam.df$variable == i] <- 2
    
  }
  
}

#Make function to bold genes which were also predicted to be CSEPs
bold_labels_cseps <- function(breaks) {
  cseps <- filter(cazyme.fam.df, variable %in% breaks) %>%
    group_by(variable) %>%
    dplyr::summarise(CSEP=mean(CSEP)) %>%
    mutate(check=str_detect(CSEP, "2")) %>% 
    pull(check)
  labels <- purrr::map2(
    breaks, cseps,
    ~ if (.y) bquote(bold(.(.x))) else bquote(plain(.(.x)))
  )
  parse(text=labels)
}

#Plot grid
gg.cazymes <- ggplot(cazyme.fam.df[!is.na(cazyme.fam.df$substrate),], aes(x=variable, y=name, fill=value)) +
  facet_grid(. ~ substrate, scales="free", space="free") +
  geom_tile(color="grey", size=0.1) +
  scale_y_discrete(limits=rev(tip.order.datedtree),
                   expand=expansion(0, 0)) +
  scale_x_discrete(labels=bold_labels_cseps) +
  coord_cartesian(clip="off",
                  ylim=c(0, 63)) +
  scale_fill_gradient(low="#F0E442", high="#CC79A7",
                      breaks=pretty_breaks(),
                      guide=guide_colourbar(title="Number of CAZymes",
                                            title.position="left",
                                            direction="horizontal",
                                            title.hjust=0,
                                            title.vjust=0.8)) +
  theme_minimal() +
  theme(legend.position="top",
        legend.direction="horizontal",
        legend.margin=margin(5, 0, -10, 0),
        legend.title=element_text(size=5, face="bold", margin=margin(0, 10, 0, 0)),
        legend.text=element_text(size=3, margin=margin(0, 3, 0, 0)),
        legend.key.size=unit(6, "pt"),
        strip.text=element_text(size=3, face="bold", margin=margin(2, 0, 2, 0)),
        panel.spacing=unit(0.15, "lines"),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=3, margin=margin(-2, 0, 0, 0)))

#Plot tree to plot beside
gg.cazymes.tree <- ggtree(dated.tree$apePhy, linetype=NA) %<+% metadata +
  geom_highlight(data=sc.df.dated, 
                 aes(node=node, fill=as.factor(box)), alpha=1, extend=300, show.legend=FALSE) +
  geom_cladelab(data=sc.df.dated,
                mapping=aes(node=node, label=sc),
                offset.text=2, offset=70, align=TRUE, barsize=0.2, fontsize=1.3) +
  coord_cartesian(clip="off") +
  scale_y_continuous(expand=expansion(0, 0),
                     limits=c(0, 63)) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tiplab(aes(label=tiplab2), fontface=tiplabel.face.datedtree, size=1.5, offset=1) +
  geom_tippoint(aes(colour=lifestyle), size=0.7, show.legend=FALSE) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(metadata$lifestyle[match(dated.tree$apePhy$tip.label,
                                                                                      metadata$name)])),
                                                 col.df$lifestyle)]) +
  new_scale_colour() +
  geom_tree(size=0.25, aes(colour=node %in% gg.datedtree$data$node[
    suppressWarnings(intersect(which(as.numeric(gg.datedtree$data$label) < 70),
                               which(gg.datedtree$data$isTip == FALSE)))]),
    show.legend=FALSE) +
  scale_colour_manual(values=c("black", "red"))

#Write to file (Supplementary Figure 7)
#tiff(file=paste0("SupplementaryFig7-", Sys.Date(), ".tiff"),
#     height=4, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(gg.cazymes.tree, gg.cazymes,
          nrow=1, rel_widths=c(1, 1), align="h", axis="bt")
#dev.off()


## REPORTED LIFESTYLES GRID ##

#Make dataframe of reported lifestyles
lifestyle.df <- metadata[c((which(colnames(metadata) == "lifestyle") + 1):ncol(metadata))]
rownames(lifestyle.df) <- metadata$name

#Order colours
col.df$lifestyle <- factor(col.df$lifestyle,
                           levels=c("plant associate", "endophyte", "plant pathogen", "saprotroph", 
                                    "animal associate", "insect mutualist",  "animal pathogen", "human pathogen", "mycoparasite"))

#Make dataframe for FFSC geographical clades
clades.df <- data.frame(node=NA,
                        clade=unique(metadata$clade.lab[metadata$clade.lab != ""]))

#Get nodes for each species complex
for (i in 1:length(clades.df$clade)) {
  clades.df$node[i] <- getMRCA(dated.tree$apePhy, metadata$name[metadata$clade.lab == clades.df$clade[i]])
}

clades.df$clade <- paste0(clades.df$clade, "\nclade")

#Plot dated tree
gg.maintree <- gg.datedtree.independent +
  new_scale_fill() +
  geom_cladelab(data=clades.df, 
                mapping=aes(node=node, label=clade),
                offset.text=2, offset=50, align=TRUE, barsize=0.2, fontsize=1.2) +
  geom_nodepoint(aes(subset=(node %in% c(64, 66))),
                 size=1, colour="black") +
  geom_text_repel(data=gg.tree.data.dated.independent[gg.tree.data.dated.independent$node %in% c(64, 66),],
                  aes(x=x-max(gg.tree.data.dated.independent$x), y=y,
                      label=c("Generic concept\nsensu lato",
                              "Generic concept\nsensu stricto")),
                  hjust=1,
                  size=1.7,
                  segment.size=0.25,
                  min.segment.length=unit(0, 'lines'),
                  xlim=c(-Inf, Inf),
                  nudge_x=-5,
                  nudge_y=2.5)

#Add lifestyle grid to dated tree
gg.lifestyles.grid <- gheatmap(gg.maintree,
                               lifestyle.df,
                               offset=85,
                               width=0.2, 
                               colnames=FALSE,
                               hjust=0,
                               color="white") +
  scale_fill_manual(values=col.df$colour[match(levels(col.df$lifestyle), col.df$lifestyle)],
                    breaks=levels(col.df$lifestyle),
                    guide=FALSE) +
  new_scale_fill() +
  geom_point(data=col.df, aes(x=1, y=1, fill=lifestyle), colour=NA, size=1, alpha=0, inherit.aes=FALSE) +
  scale_fill_manual(values=col.df$colour[match(levels(col.df$lifestyle), col.df$lifestyle)],
                    breaks=levels(col.df$lifestyle),
                    labels=str_to_sentence(levels(col.df$lifestyle)),
                    name="Lifestyle") +
  guides(fill=guide_legend(direction="horizontal",
                           title.position="left",
                           nrow=2,
                           title.hjust=0,
                           override.aes=list(shape=21, alpha=1))) +
  annotate("text", x=95, y=65, label="All reported\nlifestyles", size=2) +
  theme(plot.margin=unit(c(0, 1.7, 0, 0), "cm"),
        legend.position="top",
        legend.margin=margin(0, 0, 0, 15),
        legend.key.size=unit(5, "pt"),
        legend.text=element_text(size=5, margin=margin(0, 5, 0, 0)),
        legend.title=element_text(face="bold", size=6, margin=margin(0, 5, 0, 0)))


## NUMBER OF STRAIN SPECIFIC GENES BOXPLOT ##

#Filter presence-absence matrices for strain specific orthogroups
orthogroups.count.specific <- orthogroups.count.ingroup1[
  match(orthogroups.stats.ingroup1$orthogroup[which(orthogroups.stats.ingroup1$category == "specific")],
        rownames(orthogroups.count.ingroup1)),
  ]

CSEP.count.specific <- orthogroups.count.ingroup1[
  match(orthogroups.stats.ingroup1$orthogroup[intersect(which(orthogroups.stats.ingroup1$category == "specific"),
            which(!is.na(orthogroups.stats.ingroup1$CSEP)))],
        rownames(orthogroups.count.ingroup1)),
  ]

#Make dataframe for number of strain specific orthogroups per taxon
specific.df <- data.frame(taxon=colnames(orthogroups.count.specific),
                          CSEPs=NA,
                          Orthogroups=NA)

for (i in specific.df$taxon) {
  specific.df$CSEPs[specific.df$taxon == i] <- length(which(CSEP.count.specific[,i] > 0))
  specific.df$Orthogroups[specific.df$taxon == i] <- length(which(orthogroups.count.specific[,i] > 0))
}

#Add lifestyle
specific.df$lifestyle <- metadata$lifestyle[match(specific.df$taxon, metadata$file2)]
#Add PC1 and PC2
specific.df <- cbind(specific.df,
                     phy.pca.result.orthogroups[match(specific.df$taxon,
                                                      rownames(phy.pca.result.orthogroups)),
                                                c("PC1", "PC2")])

#Test for significant difference in number of strain specific genes between lifestyles
for (i in c("CSEPs", "Orthogroups")) {
  
  #Check for normality of residuals
  plot(ggqqplot(residuals(lm(paste(i, "~ lifestyle"), data=specific.df))) +
         ggtitle(i))
  
  #Check for homogeneity of variances
  specific.levene <- specific.df %>% levene_test(as.formula(paste(i, "~ lifestyle")))
  specific.levene <- data.frame(formula=paste0("Number of strain specific ", i, " ~ lifestyle"), specific.levene)
  assign(paste0("specific.levene.", i), specific.levene)
  
  #If homogeneity of variance assumption met...
  if (specific.levene$p > 0.05) {
    
    print("Homogeneity of variance assumption met, doing ANOVA")
    #Do ANOVA including lifestyle and first 2 PCs
    specific.anova <- specific.df %>% anova_test(as.formula(paste(i, "~ PC1 + PC2 + lifestyle")))
    #Format
    specific.anova[, c(3, 6, 7)] <- NULL
    colnames(specific.anova) <- c("Effect", "Df", "F", "p")
    specific.anova$formula <- paste0("Number of strain specific ", i, " ~ PC1 + PC2 + lifestyle")
    assign(paste0("specific.anova.", i), specific.anova)
    specific.p <- specific.anova$p[3]
    
    #If ANOVA is significant...
    if (specific.p < 0.05) {
      
      print("Number of specific ", i,
            " significant according to ANOVA, doing TukeyHSD")
      #Do TukeyHSD pairwise tests
      specific.pw.df <- specific.df %>% tukey_hsd(as.formula(paste(i, "~ lifestyle")))
      specific.pw.df[c(1, 4, 9)] <- NULL
      specific.pw.df$formula <- paste0("Number of strain specific ", i, " ~ lifestyle")
      assign(paste0("specific.pw.df.", i), specific.pw.df)
      
    } else {
      
      print(paste0("ANOVA, number of specific ", i, " not significant: p=", signif(specific.p, 1)))
      
    }
    
  } else {
    
    print("Homogeneity of variance assumption not met, doing aligned rank tranform ANOVA")
    #Do ART
    specific.anova <- aligned.rank.transform(as.formula(paste(i, "~ PC1 + PC2 + lifestyle")), data=specific.df)
    #Format
    specific.anova$significance$`Sum Sq` <- NULL
    colnames(specific.anova$significance) <- c("Df", "F", "p")
    specific.anova$significance$Effect <- rownames(specific.anova$significance)
    specific.anova$significance$formula <- paste0("Number of strain specific ", i, " ~ PC1 + PC2 + lifestyle")
    assign(paste0("specific.anova.", i), specific.anova)
    specific.p <- specific.anova$significance$p[3]
    
    #If ANOVA is significant...
    if (specific.p < 0.05) {
      
      print("Number of specific ", i,
            " significant according to aligned rank transform ANOVA, doing Games Howell test")
      #Do Games Howell pairwise tests
      specific.pw.df <- specific.df %>% games_howell_test(as.formula(paste(i, "~ lifestyle")))
      specific.pw.df[c(1, 8)] <- NULL
      specific.pw.df$formula <- paste0("Number of strain specific ", i, " ~ lifestyle")
      assign(paste0("specific.pw.df.", i), specific.pw.df)
      
    } else {
      
      print(paste0("Aligned rank transform ANOVA, number of specific ", i,
                   " not significant: p=", signif(specific.p, 1)))
      
    }
    
  }
  
}

#Make dataframe with strain sample size labels for plots
plot.labels <- as.data.frame(table(metadata$lifestyle))
colnames(plot.labels) <- c("lifestyle", "num")
plot.labels$num <- paste0("n=", plot.labels$num)

#Replace spaces with linebreaks
specific.df$lifestyle <- sub(" ", "\n", specific.df$lifestyle)
plot.labels$lifestyle <- sub(" ", "\n", plot.labels$lifestyle)

#Melt dataframe for plotting
specific.df2 <- melt(specific.df, id.vars=c("taxon", "lifestyle", "PC1", "PC2"))
#Order groups
specific.df2$variable <- factor(specific.df2$variable, levels=c("Orthogroups", "CSEPs"))

#Plot boxplot of strain specific orthogroups across lifestyles
gg.specific <- ggplot(specific.df2, aes(x=lifestyle, y=value, fill=lifestyle)) +
  facet_wrap(. ~ variable, scales="free",
             labeller=labeller(variable=c(CSEPs="CSEPs", Orthogroups="All genes"))) +
  geom_violin(fill="white", colour="grey", lty="dotted", size=0.3) +
  geom_boxplot(width=0.2, size=0.3, outlier.size=0.5) +
  geom_text(data=plot.labels,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=num),
            fontface="bold",
            hjust=0.5,
            vjust=2,
            size=1.5,
            inherit.aes=FALSE) +
  labs(y="Strain specific genes") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                      guide=FALSE) +
  scale_fill_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                    guide=FALSE) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.text=element_text(size=6, margin=margin(0, 5, 0, 0)),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        strip.text=element_text(size=6, face="bold"),
        plot.margin=margin(0, 0, 10, 0, unit="pt"))


## GENE COPY NUMBER BOXPLOT ##

#Filter presence-absence matrices for orthogroups present in at least one ingroup taxon
copynum.orthogroups <- orthogroups.count.ingroup1[which(rowSums(orthogroups.count.ingroup1) > 0),]
copynum.CSEPs <- CSEP.orthogroups
copynum.CAZymes <- CAZyme.orthogroups

#Make dataframes of gene copy number for each taxon
copynum.df.orthogroups <- data.frame(taxon=colnames(copynum.orthogroups),
                                     lifestyle=metadata$lifestyle[match(
                                       colnames(copynum.orthogroups),
                                       metadata$file2)],
                                     as.data.frame(t(copynum.orthogroups)))

copynum.df.CSEPs <- data.frame(taxon=colnames(copynum.CSEPs),
                               lifestyle=metadata$lifestyle[match(
                                 colnames(copynum.CSEPs),
                                 metadata$file2)],
                               as.data.frame(t(copynum.CSEPs)))

copynum.df.CAZymes <- data.frame(taxon=colnames(copynum.CAZymes),
                                 lifestyle=metadata$lifestyle[match(
                                   colnames(copynum.CAZymes),
                                   metadata$file2)],
                                 as.data.frame(t(copynum.CAZymes)))

#Combine dataframes for plotting
copynum.df <- rbind(data.frame(melt(copynum.df.orthogroups), gene="Orthogroups"),
                    data.frame(melt(copynum.df.CSEPs), gene="CSEPs"),
                    data.frame(melt(copynum.df.CAZymes), gene="CAZymes"))

#Remove copy number of 0
copynum.df <- copynum.df[copynum.df$value != 0,]

#Order groups
copynum.df$gene <- factor(copynum.df$gene, levels=c("Orthogroups", "CSEPs", "CAZymes"))

#Get mean copy number for each taxon
copynum.mean.df <- copynum.df  %>%
  group_by(taxon, lifestyle, gene) %>%
  summarise(mean=mean(value)) %>%
  ungroup()

#Add PC1 and PC2
copynum.mean.df <- cbind(copynum.mean.df,
                         phy.pca.result.orthogroups[match(copynum.mean.df$taxon,
                                                          rownames(phy.pca.result.orthogroups)),
                                                    c("PC1", "PC2")])

#Test for significant difference in gene copy number between lifestyles
for (i in c("CSEPs", "CAZymes", "Orthogroups")) {
  
  tmp <- copynum.mean.df[copynum.mean.df$gene == i,]
  
  #Check for normality of residuals
  plot(ggqqplot(residuals(lm(mean ~ lifestyle, data=tmp))) +
         ggtitle(i))
  
  #Check for homogeneity of variances
  copynum.levene <- tmp %>% levene_test(mean ~ lifestyle)
  copynum.levene <- data.frame(formula=paste0(i, " mean copy number ~ lifestyle"), copynum.levene)
  assign(paste0("copynum.levene.", i), copynum.levene)
  
  #If homogeneity of variance assumption met...
  if (copynum.levene$p > 0.05) {
    
    print("Homogeneity of variance assumption met, doing ANOVA")
    #Do ANOVA including lifestyle and first 2 PCs
    copynum.anova <- tmp %>% anova_test(mean ~ PC1 + PC2 + lifestyle)
    #Format
    copynum.anova[, c(3, 6, 7)] <- NULL
    colnames(copynum.anova) <- c("Effect", "Df", "F", "p")
    copynum.anova$formula <- paste0(i, " mean copy number ~ PC1 + PC2 + lifestyle")
    assign(paste0("copynum.anova.", i), copynum.anova)
    copynum.p <- copynum.anova$p[3]
    
    #If ANOVA is significant...
    if (copynum.p < 0.05) {
      
      print(i, " mean copynum significant according to ANOVA, doing TukeyHSD")
      #Do TukeyHSD pairwise tests
      copynum.pw.df <- tmp %>% tukey_hsd(mean ~ lifestyle)
      copynum.pw.df[c(1, 4, 9)] <- NULL
      copynum.pw.df$formula <- paste0(i, " mean copy number ~ + lifestyle")
      assign(paste0("copynum.pw.df.", i), copynum.pw.df)
      
    } else {
      
      print(paste0("ANOVA, ", i, " mean copynum not significant: p=", signif(copynum.p, 1)))
      
    }
    
  } else {
    
    print("Homogeneity of variance assumption not met, doing aligned rank tranform ANOVA")
    #Do ART
    copynum.anova <- aligned.rank.transform(mean ~ PC1 + PC2 + lifestyle, data=tmp)
    #Format
    copynum.anova$significance$`Sum Sq` <- NULL
    colnames(copynum.anova$significance) <- c("Df", "F", "p")
    copynum.anova$significance$Effect <- rownames(copynum.anova$significance)
    copynum.anova$significance$formula <- paste0(i, " mean copy number ~ PC1 + PC2 + lifestyle")
    assign(paste0("copynum.anova.", i), copynum.anova)
    copynum.p <- copynum.anova$significance$p[3]
    
    #If ANOVA is significant...
    if (copynum.p < 0.05) {
      
      print(i, " mean copynum significant according to aligned rank transform ANOVA, doing Games Howell test")
      #Do Games Howell pairwise tests
      copynum.pw.df <- tmp %>% games_howell_test(mean ~ lifestyle)
      copynum.pw.df[c(1, 8)] <- NULL
      copynum.pw.df$formula <- paste0(i, " mean copy number ~ lifestyle")
      assign(paste0("copynum.pw.df.", i), copynum.pw.df)
      
    } else {
      
      print(paste0("Aligned rank transform ANOVA, ", i, " mean copynum not significant: p=", signif(copynum.p, 1)))
      
    }
    
  }
  
}

#Replace spaces with linebreaks
copynum.df$lifestyle <- sub(" ", "\n", copynum.df$lifestyle)
copynum.mean.df$lifestyle <- sub(" ", "\n", copynum.mean.df$lifestyle)

#Plot boxplot of gene copy number across lifestyles
gg.copynum.mean <- ggplot(copynum.mean.df, aes(x=lifestyle, y=mean, fill=lifestyle)) +
  facet_wrap(. ~ gene,
             labeller=labeller(gene=c(CAZymes="CAZymes", CSEPs="CSEPs", Orthogroups="All genes"))) +
  geom_violin(fill="white", colour="grey", lty="dotted", size=0.3) +
  geom_boxplot(width=0.3, size=0.3, outlier.size=0.5) +
  geom_text(data=plot.labels,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=num),
            fontface="bold",
            hjust=0.5,
            vjust=2,
            size=1.5,
            show.legend=FALSE,
            inherit.aes=FALSE) +
  labs(y="Mean copy number") +
  scale_y_continuous(expand=c(0, 0),
                     limits=c(1, ceiling(max(copynum.mean.df$mean)))) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                      guide=FALSE) +
  scale_fill_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                    labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                  col.df$lifestyle)]),
                    name=NULL) +
  coord_cartesian(clip="off") +
  guides(colour=guide_legend(direction="horizontal",
                             nrow=2,
                             override.aes=list(size=1.5))) +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.text=element_text(size=6, margin=margin(0, 5, 0, 0)),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        strip.text=element_text(size=6, face="bold"),
        plot.margin=margin(5, 0, 5, 0, unit="pt"))

#Set seed for reproducible jitter
set.seed(1)

gg.copynum <- ggplot(copynum.df, aes(x=lifestyle, y=value)) +
  facet_wrap(. ~ gene,
             labeller=labeller(gene=c(CAZymes="CAZymes", CSEPs="CSEPs", Orthogroups="All genes"))) +
  geom_point(position="jitter", aes(colour=lifestyle), size=0.3) +
  geom_text(data=plot.labels,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=num),
            fontface="bold",
            hjust=0.5,
            vjust=2,
            size=1.5,
            show.legend=FALSE,
            inherit.aes=FALSE) +
  labs(y="Copy number") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                      labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                    col.df$lifestyle)]),
                      name=NULL) +
  coord_cartesian(clip="off") +
  guides(colour=guide_legend(direction="horizontal",
                             nrow=2,
                             override.aes=list(size=1.5))) +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.text=element_text(size=6, margin=margin(0, 5, 0, 0)),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        strip.text=element_text(size=6, face="bold"),
        plot.margin=margin(2, 60, 5, 50, unit="pt"))

#Get top outlier
orthogroups.stats.ingroup0[match(copynum.df$variable[which.max(copynum.df$value)], orthogroups.stats.ingroup0$orthogroup),]

#Write boxplot of strain-specific genes/copy number to file (Supplementary Figure 4)
#tiff(file=paste0("SupplementaryFig4-", Sys.Date(), ".tiff"),
#     height=4, width=4.5, unit="in", res=600, compression="lzw")
plot_grid(gg.specific,
          gg.copynum.mean,
          nrow=2,
          rel_heights=c(0.8, 1),
          labels="AUTO",
          label_size=10)
#dev.off()


## LIFESTYLE COMPARISON PERMANOVA GRID ##

#For CSEPs and all orthogroups...
for (i in c("CSEPs", "CAZymes", "orthogroups")){
  
  print(i)
  lifestyle.data <- read.csv(paste0("lifestyle_comparison/", i, "/data.csv"), row.names="genome")
  #Make distance matrix of orthogroup content
  dist <- vegdist(lifestyle.data, method="jaccard")
  #Read in lifestyle test results
  phy.pca.result <- read.csv(paste0("lifestyle_comparison/", i, "/metadata.csv"))
  #Do permanova
  permanova <- adonis2(formula=dist ~ PC1 + PC2 + lifestyle,
                       data=phy.pca.result,
                       permutations=9999)
  
  print(paste0("Phylogeny: ", round(sum(permanova$R2[1:2]) * 100), "%"))
  print(paste0("Lifestyle: ", round(sum(permanova$R2[3]) * 100), "%"))
  
  assign(paste0("permanova.", i), permanova)
  
}

#Make dataframe with pairwise PERMANOVA results
pw.lifestyle.genes <- rbind(data.frame(melt(as.matrix(read.csv(
  "lifestyle_comparison/CSEPs/pairwiseComparisons.csv",
  row.names=1)), na.rm=TRUE), data="CSEPs"),
  data.frame(melt(as.matrix(read.csv(
    "lifestyle_comparison/CAZymes/pairwiseComparisons.csv",
    row.names=1)), na.rm=TRUE), data="CAZymes"),
  data.frame(melt(as.matrix(read.csv(
    "lifestyle_comparison/orthogroups/pairwiseComparisons.csv",
    row.names=1)), na.rm=TRUE), data="Orthogroups"))

#Order groups
pw.lifestyle.genes$data <- factor(pw.lifestyle.genes$data, levels=c("Orthogroups", "CSEPs", "CAZymes"))

#Round p values
pw.lifestyle.genes$label <- round(pw.lifestyle.genes$value, digits=3)
pw.lifestyle.genes$label[which(pw.lifestyle.genes$label == 0)] <- "<0.001"

#Correct labels
pw.lifestyle.genes$Var2 <- sub("p.value.", "", pw.lifestyle.genes$Var2)
pw.lifestyle.genes$Var1 <- as.character(pw.lifestyle.genes$Var1)
pw.lifestyle.genes$Var2 <- as.character(pw.lifestyle.genes$Var2)

for (i in 1:length(pw.lifestyle.genes$Var1)) {
  pw.lifestyle.genes$Var1[i] <- as.character(col.df$lifestyle[agrep(pw.lifestyle.genes$Var1[i], col.df$lifestyle)])
  pw.lifestyle.genes$Var2[i] <- as.character(col.df$lifestyle[agrep(pw.lifestyle.genes$Var2[i], col.df$lifestyle)])
}

#Replace spaces with linebreaks
pw.lifestyle.genes$Var1 <- sub(" ", "\n", pw.lifestyle.genes$Var1)
pw.lifestyle.genes$Var2 <- sub(" ", "\n", pw.lifestyle.genes$Var2)

#Make dataframe with PERMANOVA results for labelling
permanova.df <- data.frame(lab=c(paste0("Phylogeny: ", round(sum(permanova.CSEPs$R2[1:2]) * 100),
                                        "% (p=", signif(permanova.CSEPs$`Pr(>F)`[1], 1), ")\n", 
                                        "Lifestyle: ", round(sum(permanova.CSEPs$R2[3]) * 100),
                                        "% (p=", signif(permanova.CSEPs$`Pr(>F)`[3], 1), ")"),
                                 paste0("Phylogeny: ", round(sum(permanova.CAZymes$R2[1:2]) * 100),
                                        "% (p=", signif(permanova.CAZymes$`Pr(>F)`[1], 1), ")\n", 
                                        "Lifestyle: ", round(sum(permanova.CAZymes$R2[3]) * 100),
                                        "% (p=", signif(permanova.CAZymes$`Pr(>F)`[3], 1), ")"),
                                 paste0("Phylogeny: ", round(sum(permanova.orthogroups$R2[1:2]) * 100),
                                        "% (p=", signif(permanova.orthogroups$`Pr(>F)`[1], 1), ")\n",
                                        "Lifestyle: ", round(sum(permanova.orthogroups$R2[3]) * 100), 
                                        "% (p=", signif(permanova.orthogroups$`Pr(>F)`[3], 1), ")")),
                           data=c("CSEPs", "CAZymes", "Orthogroups"))
permanova.df$data <- factor(permanova.df$data, levels=c("Orthogroups", "CSEPs", "CAZymes"))

#Plot grid of PERMANOVA p values
gg.pwperm <- ggplot(pw.lifestyle.genes, aes(Var2, Var1, fill=value>0.05)) +
  facet_grid(. ~ data, labeller=labeller(data=c(Orthogroups="All genes", CSEPs="CSEPs", CAZymes="CAZymes"))) +
  geom_tile(color="grey", size=1, alpha=0.5, show.legend=FALSE) +
  geom_text(aes(label=label, colour=value>0.05), size=1.5, show.legend=FALSE) +
  annotate("text", x=4.1, y=2, label="PERMANOVA", size=1.5, fontface="bold") +
  geom_text(data=permanova.df, aes(x=4.1, y=1.5, label=lab), size=1.5, inherit.aes=FALSE) +
  scale_fill_manual(values=c("#E69F00", "white")) +
  scale_colour_manual(values=c("black", "grey")) +
  scale_x_discrete(position="top") + 
  theme_minimal() + 
  theme(axis.text.x=element_text(colour=c("#009E73","#56B4E9", "#D55E00", "#9AE324", "dimgrey"),
                                 size=2.5, face="bold"),
        axis.text.y=element_text(colour=c("#56B4E9", "#D55E00", "#9AE324", "dimgrey", "#0072B2"),
                                 face="bold", angle=90, hjust=0.5, vjust=0, size=2.5),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text=element_text(face="bold", size=6),
        strip.placement="outside",
        panel.grid=element_blank(),
        plot.margin=unit(c(0, 0, 0, 0), "pt")) +
  coord_fixed()


## COMBINE ##

#Write to file
#tiff(file=paste0("Fig2-", Sys.Date(), ".tiff"),
#     height=8, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.lifestyles.grid, gg.gene.numbers,
                    nrow=1, rel_widths=c(6, 2), align="h", axis="bt"),
          plot_grid(gg.pwperm, gg.copynum, nrow=2, rel_heights=c(1, 1), labels=c("", "C"), label_size=10),
          nrow=2,
          rel_heights=c(2, 2),
          labels="AUTO",
          label_size=10)
dev.off()


############################    FIGURE 3    ####################################
############################   SELECTION    ####################################

## BUSTED AND aBSREL ##

#Vector of core, single-copy orthogroups
core.SC.orthogroups <- Reduce(intersect,
                              list(orthogroups.stats.ingroup0$orthogroup[which(
                                orthogroups.stats.ingroup0$copy_number == "single")],
                                orthogroups.stats.ingroup0$orthogroup[which(
                                  orthogroups.stats.ingroup0$category == "core")]))

#Core, single-copy CSEPs
core.SC.CSEPmixed <- Reduce(intersect,
                            list(orthogroups.stats.ingroup0$orthogroup[which(
                              orthogroups.stats.ingroup0$copy_number == "single")],
                              orthogroups.stats.ingroup0$orthogroup[which(
                                orthogroups.stats.ingroup0$category == "core")],
                              orthogroups.stats.ingroup0$orthogroup[which(
                                !is.na(orthogroups.stats.ingroup0$CSEP))]))

#Core, single-copy CAZymes
core.SC.CAZymemixed <- Reduce(intersect,
                              list(orthogroups.stats.ingroup0$orthogroup[which(
                                orthogroups.stats.ingroup0$copy_number == "single")],
                                orthogroups.stats.ingroup0$orthogroup[which(
                                  orthogroups.stats.ingroup0$category == "core")],
                                orthogroups.stats.ingroup0$orthogroup[which(
                                  !is.na(orthogroups.stats.ingroup0$CAZyme))]))

#Make dataframe for selection results
selection.df <- data.frame(orthogroup=core.SC.orthogroups,
                           CSEP="N", CAZyme="N", busted="N", absrel="N", busted.success="Y", absrel.success="Y")
#Orthogroups predicted to be CSEPs/CAZymes
selection.df$CSEP[match(core.SC.CSEPmixed, selection.df$orthogroup)] <- "Y"
selection.df$CAZyme[match(core.SC.CAZymemixed, selection.df$orthogroup)] <- "Y"

#Make empty list to capture aBSREL p-values for each lineage
absrel.p <- list()

#Create bar to show progress
progress.bar <- txtProgressBar(1, length(selection.df$orthogroup), initial=0, char="=", style=3)

#For each orthogroup...
for (i in selection.df$orthogroup) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, which(selection.df$orthogroup == i))
  
  #Try to read in BUSTED results
  busted.results <- tryCatch(fromJSON(paste0("selection/hyphy/busted/", i, "_BUSTED.json")), error=function(e) NULL)
  
  if (!is.null(busted.results)) {
    
    #Add whether the p-value is significant to the results dataframe
    if (busted.results$`test results`$`p-value` < 0.05)
      selection.df$busted[match(i, selection.df$orthogroup)] <- "Y"
  } else {
    selection.df$busted.success[match(i, selection.df$orthogroup)] <- "N"
  }
  
  #Try to read in aBSREL results
  absrel.results <- tryCatch(fromJSON(paste0("selection/hyphy/absrel/", i, "_aBSREL.json")), error=function(e) NULL)
  
  if (!is.null(absrel.results)) {
    
    #Add whether p-value is significant to the results list
    for (j in 1:length(names(absrel.results[["branch attributes"]][["0"]]))) {
      
      if (is.null(absrel.results[["branch attributes"]][["0"]][j][[1]][["original name"]])) {
        
        absrel.p[[names(absrel.results[["branch attributes"]][["0"]][j])]][i] <- 
          absrel.results[["branch attributes"]][["0"]][[names(absrel.results[["branch attributes"]][["0"]])[j]]][["Corrected P-value"]]
        
      } else {
        
        absrel.p[[absrel.results[["branch attributes"]][["0"]][j][[1]][["original name"]]]][i] <- 
          absrel.results[["branch attributes"]][["0"]][[names(absrel.results[["branch attributes"]][["0"]])[j]]][["Corrected P-value"]]
        
      }
      
    }
    
  } else {
    
    selection.df$absrel.success[match(i, selection.df$orthogroup)] <- "N"
    
  }
  
}

#Filter for significant aBSREL results
absrel.p.sig <- lapply(absrel.p, function(x) names(which(x < 0.05)))

#Add aBSREL results to selection dataframe 
selection.df$absrel[which(selection.df$orthogroup %in% unlist(absrel.p.sig))] <- "Y"

#Make matrices for euler plots of selection results
selection.mat.ortho <-
  c("All genes"=length(union(which(selection.df$busted.success == "Y"),
                             which(selection.df$absrel.success == "Y"))),
    "BUSTED"=0,
    "aBSREL"=0,
    "BUSTED&aBSREL"=0,
    "All genes&BUSTED"=length(Reduce(intersect,
                                     list(which(selection.df$busted == "Y"),
                                          which(selection.df$absrel == "N")))),
    "All genes&aBSREL"=length(Reduce(intersect,
                                     list(which(selection.df$busted == "N"),
                                          which(selection.df$absrel == "Y")))),
    "All genes&BUSTED&aBSREL"=length(Reduce(intersect,
                                            list(which(selection.df$busted == "Y"),
                                                 which(selection.df$absrel == "Y")))))

selection.mat.CSEPs <-
  c("CSEPs"=length(which(selection.df$CSEP == "Y")),
    "BUSTED"=0,
    "aBSREL"=0,
    "BUSTED&aBSREL"=0,
    "CSEPs&BUSTED"=length(Reduce(intersect,
                                 list(which(selection.df$busted == "Y"),
                                      which(selection.df$absrel == "N"),
                                      which(selection.df$CSEP == "Y")))),
    "CSEPs&aBSREL"=length(Reduce(intersect,
                                 list(which(selection.df$busted == "N"),
                                      which(selection.df$absrel == "Y"),
                                      which(selection.df$CSEP == "Y")))),
    "CSEPs&BUSTED&aBSREL"=length(Reduce(intersect,
                                        list(which(selection.df$busted == "Y"),
                                             which(selection.df$absrel == "Y"),
                                             which(selection.df$CSEP == "Y")))))

selection.mat.CAZymes <-
  c("CAZymes"=length(which(selection.df$CAZyme == "Y")),
    "BUSTED"=0,
    "aBSREL"=0,
    "BUSTED&aBSREL"=0,
    "CAZymes&BUSTED"=length(Reduce(intersect,
                                   list(which(selection.df$busted == "Y"),
                                        which(selection.df$absrel == "N"),
                                        which(selection.df$CAZyme == "Y")))),
    "CAZymes&aBSREL"=length(Reduce(intersect,
                                   list(which(selection.df$busted == "N"),
                                        which(selection.df$absrel == "Y"),
                                        which(selection.df$CAZyme == "Y")))),
    "CAZymes&BUSTED&aBSREL"=length(Reduce(intersect,
                                          list(which(selection.df$busted == "Y"),
                                               which(selection.df$absrel == "Y"),
                                               which(selection.df$CAZyme == "Y")))))

#Create and plot euler diagrams
for (i in c("ortho", "CSEPs", "CAZymes")) {
  
  selection.mat <- get(paste0("selection.mat.", i))
  
  if (i == "ortho") {
    
    set.seed(1)
    selection.euler <- euler(get(paste0("selection.mat.", i)), shape="ellipse")
    eulerr_options(padding=unit(0.5, "pt"))
    euler <- plot(selection.euler, 
                  fills=list(fill=c("white", "dimgrey", "dimgrey", rep("", 3), "#56B4E9"),
                             alpha=0.3),
                  labels=list(cex=0.35,
                              col="dimgrey",
                              padding=unit(c(1, 1), "pt")),
                  edges=list(col="dimgrey", lty=c("solid", "dotted", "dotted"), lwd=0.5),
                  quantities=list(col=c(rep("dimgrey", 3), "black"),
                                  cex=0.3))
    
  } else {
    
    set.seed(1)
    selection.euler <- euler(get(paste0("selection.mat.", i)), shape="ellipse")
    eulerr_options(padding=unit(0.5, "pt"))
    euler <- plot(selection.euler, 
                  fills=list(fill=c("white", "dimgrey", "dimgrey", rep("", 3), "#56B4E9"),
                             alpha=0.3),
                  labels=list(labels=c(i, rep("", 6)),
                              cex=0.35,
                              col="dimgrey",
                              padding=unit(c(1, 1), "pt")),
                  edges=list(col="dimgrey", lty=c("solid", "dotted", "dotted"), lwd=0.5),
                  quantities=list(col=c(rep("dimgrey", 3), "black"),
                                  cex=0.3))
    
  }
  
  #Convert to grob
  euler.grob <- as.grob(euler)
  
  assign(paste0("euler.grob.", i), euler.grob)
  
}

#Number of core SC orthogroups positively selected (consensus between aBSREL and BUSTED)
print(paste0("All genes: ", selection.mat.ortho[["All genes&BUSTED&aBSREL"]], ", ",
             round(selection.mat.ortho[["All genes&BUSTED&aBSREL"]] /
                     selection.mat.ortho[["All genes"]] * 100), "%"))
#Number of core SC CSEPs positively selected (consensus between aBSREL and BUSTED)
print(paste0("CSEPS: ", selection.mat.CSEPs[["CSEPs&BUSTED&aBSREL"]], ", ",
             round(selection.mat.CSEPs[["CSEPs&BUSTED&aBSREL"]] /
                     selection.mat.CSEPs[["CSEPs"]] * 100), "%"))
#Number of core SC CSEPs positively selected (consensus between aBSREL and BUSTED)
print(paste0("CAZymes: ", selection.mat.CAZymes[["CAZymes&BUSTED&aBSREL"]], ", ",
             round(selection.mat.CAZymes[["CAZymes&BUSTED&aBSREL"]] /
                     selection.mat.CAZymes[["CAZymes"]] * 100), "%"))

#Read in dataframe matching aBSREL tree nodes to RAxML-NG species tree nodes
absrel.df <- read.csv("selection/hyphy/absrel_nodes.csv")
#Add tip labels
absrel.df <- rbind(absrel.df, data.frame(node=1:length(raxmlng$tip.label), absrel=raxmlng$tip.label))

##Add number of significant aBSREL tests per taxon that are also supported by busted
absrel.p.count <- lengths(lapply(absrel.p.sig, function(x)
  x[which(x %in% selection.df$orthogroup[which(selection.df$busted == "Y")])]))
absrel.df$num <- absrel.p.count[match(absrel.df$absrel, names(absrel.p.count))]

#Add columns for CSEP and CAZyme labels
absrel.df[c("CSEP", "CAZyme", "CAZyme_family", "CSEP_name")] <- NA

#For CSEPs and CAZymes...
for (i in c("CSEP", "CAZyme")) {
  
  #Add IDs for positively selected genes
  absrel.p.names <- unlist(lapply(absrel.p.sig, function(x)
    paste(x[intersect(which(x %in% get(paste0("core.SC.", i, "mixed"))),
                      which(x %in% selection.df$orthogroup[which(selection.df$busted == "Y")]))],
          collapse=" ")))
  absrel.p.names <- absrel.p.names[absrel.p.names != ""]
  absrel.df[,i] <- absrel.p.names[match(absrel.df$absrel, names(absrel.p.names))]
  
  if (i == "CSEP") {
    
    suffix <- "_name"
    
  } else if (i == "CAZyme") {
    
    suffix <- "_family"
    
  }
  
  #Add gene names if possible
  for (j in 1:length(absrel.df[,i])) {
    
    tmp <- unlist(str_split(absrel.df[j, i], " "))
    
    if (length(tmp) > 1) {
      
      absrel.df[j, paste0(i, suffix)] <-
        paste(unique(na.omit(orthogroups.stats.ingroup0[
          match(tmp, orthogroups.stats.ingroup0$orthogroup),
          paste0(i, suffix)
        ])),
        collapse=" ")
      
    } else {
      
      absrel.df[j, paste0(i, suffix)] <- 
        orthogroups.stats.ingroup0[match(absrel.df[j, i],
                                         orthogroups.stats.ingroup0$orthogroup),
                                   paste0(i, suffix)]
      
    }
    
  }
  
  absrel.df[,i] <- gsub("OG000", "", absrel.df[,i])
  absrel.df[,i][absrel.df[,i] == "<NA>"] <- NA
  absrel.df[,paste0(i, suffix)][absrel.df[,paste0(i, suffix)] == ""] <- NA
  
  #Wrap gene IDs and names
  absrel.df[,i]  <- gsub('(?=(?:.{10})+$)', "\n", absrel.df[,i], perl=TRUE)
  
  #for (k in 1:length(absrel.df[!is.na(absrel.df[paste0(i, suffix)]), paste0(i, suffix)])) {
  #  
  #  absrel.df[!is.na(absrel.df[paste0(i, suffix)]), paste0(i, suffix)][k] <- 
  #    paste0(stri_wrap(absrel.df[!is.na(absrel.df[paste0(i, suffix)]), paste0(i, suffix)][k],
  #                     nchar(absrel.df[!is.na(absrel.df[paste0(i, suffix)]), paste0(i, suffix)][k])/2,
  #                     whitespace_only=TRUE), collapse="\n")
  #  
  #}
  
  assign(paste0("absrel.p.", i), absrel.p.names)
  
}

#Plot species tree with no branch lengths
gg.selection <- ggtree(raxmlng.tree, branch.length="none", linetype=NA) %<+% metadata

#Create vector with the order of the tree tips for plotting
tip.order.speciestree <- rev(gg.tree.data.raxmlng$label[gg.tree.data.raxmlng$isTip == "TRUE"])
#Create vector of fontface for new genomes in this study
label.face.speciestree <- rev(metadata$new[match(tip.order.speciestree, metadata$name)])
for (i in 1:length(label.face.speciestree)){
  if (label.face.speciestree[i] == "Y") {
    label.face.speciestree[i] <- "bold.italic"
  } else {
    label.face.speciestree[i] <- "italic"
  }
}
#Fontface vector
tiplabel.face.speciestree <- rep("italic", length(raxmlng.tree$tip.label))
tiplabel.face.speciestree[which(metadata$new[match(raxmlng.tree$tip.label, metadata$name)] == "Y")] <- "bold.italic"

#Make dataframe of species complex nodes
sc.df.spec <- data.frame(sc=unique(metadata$speciescomplex.abb),
                         node=NA)

#Get nodes for each species complex
for (i in 1:length(sc.df.spec$sc)) {
  sc.df.spec$node[i] <- MRCA(raxmlng.tree, metadata$name[metadata$speciescomplex.abb == sc.df.spec$sc[i]])
}

#Make alternated coding for highlights on tree
sc.df.spec <- sc.df.spec[match(na.omit(unique(gg.tree.data.raxmlng$speciescomplex.abb)), sc.df.spec$sc),]
sc.df.spec$box <- rep(c(0,1), length.out=length(sc.df.spec$sc))

#Add species complex highlights
gg.selection <- gg.selection +
  geom_highlight(data=sc.df.spec, 
                 aes(node=node, fill=as.factor(box)), alpha=1, extend=10, show.legend=FALSE) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC"))

#Write numbered selection tree to file (Supplementary Figure 5)
#tiff(file=paste0("SupplementaryFig5-", Sys.Date(), ".tiff"),
#     height=4, width=6.75, unit="in", res=600, compression="lzw")

gg.selection %<+% absrel.df +
  geom_tree(size=1, aes(colour=num)) +
  xlim(0, 21) +
  scale_colour_gradient2(trans="pseudo_log",
                         low="black", mid="#F0E442", high="#CC79A7", midpoint=1,
                         breaks=c(0, 10, 100),
                         limits=c(0, 100),
                         labels=c(0 , 10, 100),
                         na.value="dimgrey",
                         guide=guide_colourbar(title="Genes showing\npositive selection",
                                               title.position="top",
                                               direction="horizontal")) +
  new_scale_fill() +
  geom_label2(aes(x=branch, label=num, fill=num),
              size=1, label.size=0, label.padding=unit(1.5, "pt"), fontface="plain") +
  scale_fill_gradient2(trans="pseudo_log",
                       low="black", mid="#F0E442", high="#CC79A7", midpoint=1,
                       breaks=c(0, 10, 100),
                       limits=c(0, 100),
                       labels=c(0 , 10, 100),
                       na.value="dimgrey",
                       guide=guide_colourbar(title="Genes showing\npositive selection",
                                             title.position="top",
                                             direction="horizontal")) +
  geom_tiplab(aes(label=tiplab2), fontface=tiplabel.face.speciestree, size=1.5, offset=0.05, colour="black") +
  coord_cartesian(clip="off") +
  annotation_custom(euler.grob.ortho, xmin=0, xmax=6, ymin=45, ymax=62) +
  annotation_custom(euler.grob.CSEPs, xmin=3, xmax=10, ymin=50, ymax=62) +
  annotation_custom(euler.grob.CAZymes, xmin=3, xmax=8, ymin=42, ymax=50) +
  theme(legend.box="vertical",
        legend.title=element_text(size=5, face="bold"),
        legend.text=element_text(size=5, margin=margin(0, 3, 0, 0)),
        legend.position=c(0.13, 0.45),
        legend.key.size=unit(6, "pt"),
        legend.spacing=unit(0, "pt"),
        plot.margin=margin(0, 0, 5, -16, unit="pt"))

#dev.off()

#Make grob for gene label legend
genelabels.grob <- as.grob(ggplot() +
                             geom_label(aes(x=1, y=2),
                                        label="CSEP",
                                        fill="dimgrey",
                                        colour="white",
                                        size=1.5, label.size=0,
                                        label.padding=unit(2.25, "pt"), fontface="bold") +
                             geom_label(aes(x=1, y=1),
                                        label="CAZyme",
                                        fill="dimgrey",
                                        colour="white",
                                        size=1.5, label.size=0,
                                        label.padding=unit(2.25, "pt"), fontface="italic") +
                             coord_cartesian(clip="off") +
                             theme_void())

#Plot aBSREL results on tree
gg.selection2 <- gg.selection %<+% absrel.df +
  geom_tree(size=1, aes(colour=num)) +
  xlim(0, 21) +
  scale_colour_gradient2(trans="pseudo_log",
                         low="black", mid="#F0E442", high="#CC79A7", midpoint=1,
                         breaks=c(0, 10, 100),
                         limits=c(0, 100),
                         labels=c(0 , 10, 100),
                         na.value="dimgrey",
                         guide=guide_colourbar(title="Genes showing\npositive selection",
                                               title.position="top",
                                               direction="horizontal")) +
  geom_label_repel(aes(x=ifelse(!is.na(CAZyme), branch-0.4, branch), 
                       label=ifelse(branch < 16, CSEP_name, "")),
                   nudge_x=0, nudge_y=0.5,
                   colour="dimgrey",
                   size=1,
                   min.segment.length=unit(0, 'lines'),
                   label.padding=unit(1.5, "pt"), fontface="bold") +
  geom_label_repel(aes(x=ifelse(!is.na(CAZyme), branch-0.4, branch), 
                       label=ifelse(branch > 16, CSEP_name, "")),
                   nudge_x=-0.05, nudge_y=0.5,
                   colour="dimgrey",
                   size=1,
                   min.segment.length=unit(0, 'lines'),
                   label.padding=unit(1.5, "pt"), fontface="bold") +
  geom_label2(aes(x=ifelse(!is.na(CAZyme), branch-0.4, branch), label=CSEP, fill=num),
              fill="dimgrey",
              colour="white",
              size=1, 
              label.size=0,
              label.padding=unit(1.5, "pt"),
              fontface="bold") +
  geom_label_repel(aes(x=ifelse(!is.na(CSEP), branch+0.4, branch),
                       label=CAZyme_family),
                   nudge_x=-0, nudge_y=0.5,
                   colour="dimgrey",
                   size=1,
                   min.segment.length=unit(0, 'lines'),
                   label.padding=unit(1.5, "pt"), fontface="italic") +
  #geom_label_repel(aes(x=ifelse(!is.na(CSEP), branch+0.4, branch),
  #                    label=ifelse(branch > 17, CAZyme_family, "")),
  #                 nudge_x=-3,
  #                 direction="x",
  #                 colour="dimgrey",
  #                 size=1,
  #                 min.segment.length=unit(0, 'lines'),
  #                 label.padding=unit(1.5, "pt"), fontface="italic") +
  geom_label2(aes(x=ifelse(!is.na(CSEP), branch+0.4, branch), label=CAZyme, colour=num),
              fill="dimgrey",
              colour="white",
              size=1, 
              label.size=0,
              label.padding=unit(1.5, "pt"),
              fontface="italic") +
  geom_tiplab(aes(label=tiplab2), fontface=tiplabel.face.speciestree, size=1.5, offset=0.05, colour="black") +
  new_scale_colour() +
  geom_tippoint(aes(colour=lifestyle), size=0.8) +
  scale_colour_manual(values=col.df$colour[
    na.omit(match(sort(unique(metadata$lifestyle[match(iqtree.tree$tip.label, metadata$name)])),
                  col.df$lifestyle))],
                      labels=str_to_sentence(col.df$lifestyle[
                        na.omit(match(sort(unique(metadata$lifestyle[match(iqtree.tree$tip.label, metadata$name)])),
                                      col.df$lifestyle))]),
                      na.translate=FALSE,
                      guide=guide_legend(title="Lifestyle",
                                         ncol=2,
                                         title.hjust=0,
                                         order=1)) +
  coord_cartesian(clip="off") +
  annotation_custom(euler.grob.ortho, xmin=-0.5, xmax=5.5, ymin=43, ymax=60) +
  annotation_custom(euler.grob.CSEPs, xmin=2.5, xmax=9.5, ymin=48, ymax=60) +
  annotation_custom(euler.grob.CAZymes, xmin=2.5, xmax=7.5, ymin=40, ymax=48) +
  annotation_custom(genelabels.grob, xmin=0.1, xmax=1.1, ymin=19, ymax=22) +
  theme(legend.box="vertical",
        legend.margin=margin(10, 0, 0, 0),
        legend.title=element_text(size=5, face="bold"),
        legend.text=element_text(size=5, margin=margin(0, 3, 0, 0)),
        legend.position=c(0.13, 0.52),
        legend.key.size=unit(6, "pt"),
        legend.spacing=unit(0, "pt"),
        plot.margin=margin(0, 0, 5, -16, unit="pt"))

#Number of CSEPs positively selected on external branches
length(which(!is.na(gg.selection2$data$CSEP[gg.selection2$data$isTip == TRUE])))
#Number of CAZymes positively selected on external branches
length(which(!is.na(gg.selection2$data$CAZyme[gg.selection2$data$isTip == TRUE])))

#Test for significant difference in number of positively selected genes on external branches between lifestyles
selection.tmp <- gg.selection2$data[which(gg.selection2$data$num[gg.selection2$data$isTip == TRUE] > 0),
                                    c("short.tip", "lifestyle", "num")]
#Add PC 1 and 2
selection.tmp <- cbind(selection.tmp,
                       phy.pca.result.orthogroups[match(selection.tmp$short.tip, phy.pca.result.orthogroups$genome), c("PC1", "PC2")])

#Check for normality of residuals
plot(ggqqplot(residuals(lm(num ~ lifestyle, data=selection.tmp))) +
       ggtitle(i))

#Check for homogeneity of variances
selection.levene <- selection.tmp %>% levene_test(num ~ lifestyle)
selection.levene <- data.frame(formula="num ~ lifestyle", selection.levene)

#If homogeneity of variance assumption met...
if (selection.levene$p > 0.05) {
  
  print("Homogeneity of variance assumption met, doing ANOVA")
  #Do ANOVA including lifestyle and first 2 PCs
  selection.anova <- selection.tmp %>% anova_test(num ~ PC1 + PC2 + lifestyle)
  #Format
  selection.anova[, c(3, 6, 7)] <- NULL
  selection.anova$ges <- NULL
  colnames(selection.anova) <- c("Effect", "Df", "F", "p")
  selection.anova$formula <- "num ~ lifestyle"
  
  #If ANOVA is significant...
  if (selection.anova$p[3] < 0.05) {
    
    print(paste0("Number of positively selected genes significant according to ANOVA, doing TukeyHSD"))
    #Do TukeyHSD pairwise tests
    selection.pw.df <- selection.tmp %>% tukey_hsd(num ~ lifestyle)
    selection.pw.df[c(1, 4, 9)] <- NULL
    selection.pw.df$formula <- "num ~ lifestyle"
    
  } else {
    
    print(paste0("ANOVA, number of positively selected genes not significant: p=", signif(selection.anova$p[3], 1)))
    
  }
  
} else {
  
  print("Homogeneity of variance assumption not met, doing aligned rank tranform ANOVA")
  #Do ART
  selection.anova <- aligned.rank.transform(num ~ lifestyle, data=selection.tmp)
  #Format
  selection.anova$significance$`Sum Sq` <- NULL
  colnames(selection.anova$significance) <- c("Df", "F", "p")
  selection.anova$significance$Effect <- rownames(selection.anova$significance)
  selection.anova$significance$formula <- "num ~ lifestyle"
  
  #If ANOVA is significant...
  if (selection.anova$significance$p[3] < 0.05) {
    
    print(paste0("Number of positively selected genes significant according to aligned rank transform ANOVA,
                 doing Games Howell test"))
    #Do Games Howell pairwise tests
    selection.pw.df <- selection.tmp %>% games_howell_test(num ~ lifestyle)
    selection.pw.df[c(1, 8)] <- NULL
    selection.pw.df$formula <- "num ~ lifestyle"
    
  } else {
    
    print(paste0("Aligned rank transform ANOVA, number of positively selected genes not significant: p=",
                 signif(selection.anova$significance$p[3], 1)))
    
  }
  
}


## CONTRAST-FEL VIOLIN PLOTS ##

#Make dataframe for Contrast-FEL results
contrastfel.df <- data.frame(ortho=rep(core.SC.orthogroups, each=length(na.omit(unique(metadata$lifestyle)))),
                             lifestyle=rep(na.omit(unique(metadata$lifestyle)), times=length(core.SC.orthogroups)),
                             increase=NA,
                             decrease=NA)

contrastfel.df$CSEP <- orthogroups.stats.ingroup1$CSEP[match(contrastfel.df$ortho,
                                                             orthogroups.stats.ingroup1$orthogroup)]
contrastfel.df$CAZyme <- orthogroups.stats.ingroup1$CAZyme[match(contrastfel.df$ortho,
                                                                 orthogroups.stats.ingroup1$orthogroup)]

#Create bar to show progress
progress.bar <- txtProgressBar(1, length(contrastfel.df$ortho), initial=0, char="=", style=3)

counter <- 0

#For each orthogroup...
for (i in core.SC.orthogroups) {
  
  #For each lifestyle...
  for (j in na.omit(unique(metadata$lifestyle))) {
    
    counter <- counter + 1
    
    #Update progress bar
    setTxtProgressBar(progress.bar, counter)
    
    lifestyle <- sub(" ", "", j)
    
    #Try to read in BUSTED results
    contrastfel.results <- tryCatch(fromJSON(paste0("selection/hyphy/contrast-fel/",
                                                    i, "_Contrast-FEL_", lifestyle, ".json")),
                                    error=function(e) NULL)
    
    if (!is.null(contrastfel.results)) {
      
      #Get site rates for background branches
      background.rate <- contrastfel.results$MLE$content$`0`[,2][
        which(contrastfel.results$MLE$content$`0`[,6] <= 0.2)] 
      #Get site rates for test-set branches
      lifestyle.rate <- contrastfel.results$MLE$content$`0`[,3][
        which(contrastfel.results$MLE$content$`0`[,6] <= 0.2)] 
      
      #Number of sites with higher rates than background branches
      contrastfel.df$increase[intersect(which(contrastfel.df$ortho == i), which(contrastfel.df$lifestyle == j))] <- 
        length(which(background.rate < lifestyle.rate))
      #Number of sites with lower rates than background branches
      contrastfel.df$decrease[intersect(which(contrastfel.df$ortho == i), which(contrastfel.df$lifestyle == j))] <- 
        length(which(background.rate > lifestyle.rate))
      
    }
    
  }
  
}

#Test for significant difference in number of sites with different relative rate between lifestyles
for (i in c("increase", "decrease")) {
  
  tmp <- contrastfel.df[which(contrastfel.df[,i] > 0),]
  
  #Check for normality of residuals
  plot(ggqqplot(residuals(lm(paste(i, "~ lifestyle"), data=tmp))) +
         ggtitle(i))
  
  #Check for homogeneity of variances
  contrastfel.levene <- tmp %>% levene_test(as.formula(paste(i, "~ lifestyle")))
  contrastfel.levene <- data.frame(formula=paste0(i, " ~ lifestyle"), contrastfel.levene)
  assign(paste0("contrastfel.levene.", i), contrastfel.levene)
  
  #If homogeneity of variance assumption met...
  if (contrastfel.levene$p > 0.05) {
    
    print("Homogeneity of variance assumption met, doing ANOVA")
    #Do ANOVA including lifestyle and first 2 PCs
    contrastfel.anova <- tmp %>% anova_test(as.formula(paste(i, "~ lifestyle")))
    #Format
    contrastfel.anova[, c(3, 6, 7)] <- NULL
    contrastfel.anova$ges <- NULL
    colnames(contrastfel.anova) <- c("Effect", "Df", "F", "p")
    contrastfel.anova$formula <- paste0(i, " ~ lifestyle")
    assign(paste0("contrastfel.anova.", i), contrastfel.anova)
    contrastfel.p <- contrastfel.anova$p
    
    #If ANOVA is significant...
    if (contrastfel.p < 0.05) {
      
      print(paste0("Number of sites with ", i,
                   "d sites significant according to ANOVA, doing TukeyHSD"))
      #Do TukeyHSD pairwise tests
      contrastfel.pw.df <- tmp %>% tukey_hsd(as.formula(paste(i, "~ lifestyle")))
      contrastfel.pw.df[c(1, 4, 9)] <- NULL
      contrastfel.pw.df$formula <- paste0(i, " ~ lifestyle")
      assign(paste0("contrastfel.pw.df.", i), contrastfel.pw.df)
      
    } else {
      
      print(paste0("ANOVA, ", i, "d sites not significant: p=", signif(contrastfel.p, 1)))
      
    }
    
  } else {
    
    print("Homogeneity of variance assumption not met, doing aligned rank tranform ANOVA")
    #Do ART
    contrastfel.anova <- aligned.rank.transform(as.formula(paste(i, "~ lifestyle")), data=tmp)
    #Format
    contrastfel.anova$significance$`Sum Sq` <- NULL
    colnames(contrastfel.anova$significance) <- c("Df", "F", "p")
    contrastfel.anova$significance$Effect <- rownames(contrastfel.anova$significance)
    contrastfel.anova$significance$formula <- paste0(i, " ~ lifestyle")
    assign(paste0("contrastfel.anova.", i), contrastfel.anova)
    contrastfel.p <- contrastfel.anova$significance$p
    
    #If ANOVA is significant...
    if (contrastfel.p < 0.05) {
      
      print(paste0("Number of sites with ", i,
                   "d sites significant according to aligned rank transform ANOVA, doing Games Howell test"))
      #Do Games Howell pairwise tests
      contrastfel.pw.df <- tmp %>% games_howell_test(as.formula(paste(i, "~ lifestyle")))
      contrastfel.pw.df[c(1, 8)] <- NULL
      contrastfel.pw.df$formula <- paste0(i, " ~ lifestyle")
      assign(paste0("contrastfel.pw.df.", i), contrastfel.pw.df)
      
    } else {
      
      print(paste0("Aligned rank transform ANOVA, ", i, "d sites not significant: p=", signif(contrastfel.p, 1)))
      
    }
    
  }
  
}

#Make dataframe for labelling plot
sitelabels.df <-
  rbind(data.frame(test=
                     multcompLetters(setNames(contrastfel.pw.df.increase$p.adj,
                                              paste0(contrastfel.pw.df.increase$group1,
                                                     "-", contrastfel.pw.df.increase$group2)))$Letters,
                   lifestyle=
                     names(multcompLetters(setNames(contrastfel.pw.df.increase$p.adj,
                                                    paste0(contrastfel.pw.df.increase$group1,
                                                           "-", contrastfel.pw.df.increase$group2)))$Letters),
                   variable="increase"),
        data.frame(test=NA,
                   lifestyle=names(multcompLetters(setNames(contrastfel.pw.df.increase$p.adj,
                                                            paste0(contrastfel.pw.df.increase$group1,
                                                                   "-", contrastfel.pw.df.increase$group2)))$Letters),
                   variable="decrease"))

#Add sample size
sitelabels.df$num <- NA

for (i in sitelabels.df$lifestyle) {
  
  sitelabels.df$num[intersect(which(sitelabels.df$lifestyle == i),
                              which(sitelabels.df$variable == "increase"))] <- 
    length(intersect(which(contrastfel.df$lifestyle == i),
                     which(contrastfel.df$increase > 0)))
  
  sitelabels.df$num[intersect(which(sitelabels.df$lifestyle == i),
                              which(sitelabels.df$variable == "decrease"))] <- 
    length(intersect(which(contrastfel.df$lifestyle == i),
                     which(contrastfel.df$decrease > 0)))
  
}

#Replace spaces with linebreaks
sitelabels.df$lifestyle <- sub(" ", "\n", sitelabels.df$lifestyle)

#Number of genes with a different pressure relative to other lifestyles
length(unique(contrastfel.df$ortho[union(which(contrastfel.df$increase > 0), which(contrastfel.df$decrease > 0))]))

#Genes with both higher and lower relatives rates
contrastfel.df[intersect(which(contrastfel.df$increase > 0), which(contrastfel.df$decrease > 0)),]

#Melt dataframe for plotting
contrastfel.sites <- melt(contrastfel.df)
#Remove genes with no sites
contrastfel.sites <- contrastfel.sites[which(contrastfel.sites$value > 0),]

#Add column for highlighting positively selected CSEPs/CAZymes
contrastfel.sites$positive.selection <- NA

contrastfel.sites$positive.selection[
  contrastfel.sites$ortho %in% unique(contrastfel.sites$ortho[
    union(which(!is.na(contrastfel.sites$CSEP)),
          which(!is.na(contrastfel.sites$CAZyme)))
  ][
    contrastfel.sites$ortho[
      union(which(!is.na(contrastfel.sites$CSEP)),
            which(!is.na(contrastfel.sites$CAZyme)))
    ] %in% selection.df$orthogroup[
      intersect(which(selection.df$busted == "Y"),
                which(selection.df$absrel == "Y"))
    ]
  ])
] <- "Y"


#Remove prefix
contrastfel.sites$ortho <- sub("OG000", "", contrastfel.sites$ortho)
#Replace spaces with linebreaks
contrastfel.sites$lifestyle <- sub(" ", "\n", contrastfel.sites$lifestyle)
#Subset data for labels
contrastfel.sites.labels <- contrastfel.sites[which(contrastfel.sites$positive.selection == "Y"),]

#Plot violin plot of Contrast-FEL results
gg.siterates <- ggplot(contrastfel.sites, aes(x=lifestyle, y=value, fill=lifestyle)) +
  facet_wrap(. ~ variable, labeller=labeller(variable=c(increase="Higher relative selective pressure",
                                                        decrease="Lower relative selective pressure"))) +
  geom_violin(aes(colour=lifestyle), size=0.3) +
  geom_point(data=contrastfel.sites.labels,
             size=0.2,
             colour="black") +
  geom_text(data=sitelabels.df,
            aes(x=lifestyle, y=Inf, label=test),
            family="mono",
            hjust=0.5,
            vjust=2,
            size=2,
            inherit.aes=FALSE) +
  geom_text(data=sitelabels.df,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=paste0("n=", num)),
            fontface="bold",
            hjust=0.5,
            vjust=6,
            size=1.5,
            inherit.aes=FALSE) +
  geom_label_repel(data=contrastfel.sites.labels,
                   aes(label=ortho,
                       fontface=ifelse(!is.na(CSEP), "bold", "italic")),
                   fill="dimgrey",
                   colour="white",
                   segment.size=0.2,
                   segment.colour="grey",
                   min.segment.length=unit(0, 'lines'),
                   label.padding=unit(0.1, "lines"),
                   nudge_x=0.5,
                   direction="y",
                   max.overlaps=Inf,
                   size=1.5) +
  labs(y="Number of sites") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.15)),
                     limits=c(0, max(contrastfel.sites$value))) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                      labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                    col.df$lifestyle)]),
                      name=NULL) +
  scale_fill_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                    labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                  col.df$lifestyle)]),
                    name=NULL) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(legend.position="none",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(6, "pt"),
        legend.text=element_text(size=4, margin=margin(0, 5, 0, 0)),
        axis.text.x=element_text(colour=c("#009E73","#56B4E9", "#D55E00", "#9AE324", "dimgrey", "#0072B2"),
                                 size=4, face="bold", vjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        strip.text=element_text(size=6, face="bold"),
        plot.margin=margin(10, 10, 10, 10, unit="pt"))

## COMBINE ##

#Write to file
#tiff(file=paste0("Fig3-", Sys.Date(), ".tiff"),
#     height=6, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(gg.selection2,
          gg.siterates,
          nrow=2,
          rel_heights=c(2, 1),
          labels="AUTO",
          label_size=10)
#dev.off()


##########################    FIGURES 4 & 5    #################################
########################   CODON OPTIMISATION    ###############################

#Load codon optimisation results
load("selection/codon_optimisation/codon_optimisation-2022-02-14.RData")

## LIFESTYLE CODON OPTIMISATION BOXPLOT ##

#Add lifestyle, name and PCs to optimisation (S value) dataframe
s.df$lifestyle <- metadata$lifestyle[match(s.df$taxon, metadata$file2)]
s.df$lifestyle <- sub(" ", "\n", s.df$lifestyle)
s.df$name <- metadata$name[match(s.df$taxon, metadata$file2)]
s.df <- cbind(s.df,
              phy.pca.result.orthogroups[match(s.df$taxon,
                                               rownames(phy.pca.result.orthogroups)), c("PC1", "PC2")])
s.df <- s.df %>% select(name, everything())

#Remove outgroup
s.df <- s.df[-which(is.na(s.df$lifestyle)),]

#Test for significant difference in S values between lifestyles
#Check for normality of residuals
plot(ggqqplot(residuals(lm(S ~ lifestyle, data=s.df))))

#Check for homogeneity of variances
slifestyles.levene <- s.df %>% levene_test(S ~ lifestyle)
slifestyles.levene <- data.frame(formula="S ~ lifestyle", slifestyles.levene)

#If homogeneity of variance assumption met...
if (slifestyles.levene$p > 0.05) {
  
  print("Homogeneity of variance assumption met, doing ANOVA")
  #Do ANOVA
  slifestyles.anova <- s.df %>% anova_test(S ~ PC1 + PC2 + lifestyle)
  #Format
  slifestyles.anova[, c(3, 6, 7)] <- NULL
  colnames(slifestyles.anova) <- c("Effect", "Df", "F", "p")
  slifestyles.anova$formula <- "S ~ PC1 + PC2 + lifestyle"
  slifestyle.p <- slifestyles.anova$p[3]
  
  #If ANOVA is significant...
  if (slifestyle.p < 0.05) {
    
    print("Lifestyle significant according to ANOVA, doing TukeyHSD")
    #Do TukeyHSD pairwise tests
    slifestyles.pw.df <- s.df %>% tukey_hsd(S ~ lifestyle)
    slifestyles.pw.df$formula <- "S ~ lifestyle"
    
  } else {
    
    print(paste0("ANOVA lifestyle not significant: p=", slifestyle.p))
    
  }
  
} else {
  
  print("Homogeneity of variance assumption not met, doing aligned rank tranform ANOVA")
  #Do ART
  slifestyles.anova <- aligned.rank.transform(S ~ PC1 + PC2 + lifestyle, data=s.df)
  #Format
  slifestyles.anova$significance$`Sum Sq` <- NULL
  colnames(slifestyles.anova$significance) <- c("Df", "F", "p")
  slifestyles.anova$significance$Effect <- rownames(slifestyles.anova$significance)
  slifestyles.anova$significance$formula <- "S ~ PC1 + PC2 + lifestyle"
  assign(paste0("slifestyles.anova.", i), slifestyles.anova)
  slifestyle.p <- slifestyles.anova$significance$p[3]
  
  #If ANOVA is significant...
  if (slifestyle.p < 0.05) {
    
    print("Lifestyle significant according to aligned rank transform ANOVA, doing Games Howell test")
    #Do Games Howell pairwise tests
    slifestyles.pw.df <- s.df %>% games_howell_test(S ~ lifestyle)
    slifestyles.pw.df$formula <- "S ~ lifestyle"
    
  } else {
    
    print(paste0("Aligned rank transform ANOVA, lifestyle not significant: p=", signif(slifestyle.p, 1)))
    
  }
  
}

#Make dataframe for labelling plot
codonlabelslifestyles.df <-
  data.frame(test=
               multcompLetters(setNames(slifestyles.pw.df$p.adj,
                                        paste0(slifestyles.pw.df$group1,
                                               "-", slifestyles.pw.df$group2)))$Letters,
             lifestyle=
               names(multcompLetters(setNames(slifestyles.pw.df$p.adj,
                                              paste0(slifestyles.pw.df$group1,
                                                     "-", slifestyles.pw.df$group2)))$Letters))

#Melt dataframe
s.df2 <- melt(s.df[-which(colnames(s.df) == "S")],
              id.vars=c("name", "taxon", "lifestyle", "PC1", "PC2"),
              variable.name="gene_type")

levene.list <- list()
anova.list <- list()

#Test for significant difference in S values between gene types for each lifestyle
for (i in unique(s.df$lifestyle)) {
  
  print(i)
  tmp <- s.df2[s.df2$lifestyle == i,]
  
  #If group is at least 6...
  if (length(unique(tmp$name)) > 3) {
    
    #Check for homogeneity of variances
    sgenes.levene <- tmp %>% levene_test(value ~ gene_type)
    sgenes.levene <- data.frame(formula=paste0(i, ": S ~ gene type"), sgenes.levene)
    levene.list[i] <- sgenes.levene$p
    assign(paste0("sgenes.levene.", gsub("\n", "", i)), sgenes.levene)
    
    #If homogeneity of variance assumption met...
    if (sgenes.levene$p > 0.05) {
      
      print("Homogeneity of variance assumption met, doing ANOVA")
      #Do ANOVA
      sgenes.anova <- tmp %>% anova_test(value ~ PC1 + PC2 + gene_type)
      #Format
      sgenes.anova[, c(3, 6, 7)] <- NULL
      colnames(sgenes.anova) <- c("Effect", "Df", "F", "p")
      sgenes.anova$formula <- paste0(i, ": S ~ PC1 + PC2 + gene type")
      assign(paste0("sgenes.anova.", gsub("\n", "", i)), sgenes.anova)
      anova.list[i] <- sgenes.anova$p[3]
      
    } else {
      
      print("Homogeneity of variance assumption not met, doing aligned rank tranform ANOVA")
      #Do ART
      sgenes.anova <- aligned.rank.transform(value ~ PC1 + PC2 + gene_type, data=tmp)
      #Format
      sgenes.anova$significance$`Sum Sq` <- NULL
      colnames(sgenes.anova$significance) <- c("Df", "F", "p")
      sgenes.anova$significance$Effect <- rownames(sgenes.anova$significance)
      sgenes.anova$significance$formula <- paste0(i, ": S ~ PC1 + PC2 + gene type")
      assign(paste0("sgenes.anova.", gsub("\n", "", i)), sgenes.anova)
      anova.list[i] <- sgenes.anova$significance$p[3]
      
    } 
    
  } else {
    
    print("Too few data points for testing")
    
  }
  
}

#Remove non-significant lifestyles/those with too few data points for pairwise comparison
s.df3 <- s.df2[s.df2$lifestyle %in% names(anova.list),]

if (TRUE %in% (anova.list > 0.05)) {
  
  s.df3 <- s.df3[s.df3$lifestyle %in% names(anova.list)[which(anova.list < 0.05)],]
  
}

#If all ANOVAs were significant...
if (FALSE %in% (levene.list > 0.05)) {
  
  print("Some significant lifestyles did not have homogeneity of variance, doing Games Howell tests")
  #Do Games Howell pairwise tests
  sgenes.pw.df <- s.df3 %>%
    group_by(lifestyle) %>%
    games_howell_test(value ~ gene_type)
  sgenes.pw.df$formula <- paste0(sgenes.pw.df$lifestyle, ": S ~ gene type")
  
} else {
  
  print("All significant lifestyle had homogeneity of variance, doing TukeyHSD")
  
  #Do TukeyHSD pairwise tests
  sgenes.pw.df <-  s.df3 %>%
    group_by(lifestyle) %>%
    tukey_hsd(value ~ gene_type)
  sgenes.pw.df$formula <- paste0(sgenes.pw.df$lifestyle, ": S ~ gene type")
  
}

#Plot boxplots of codon optimisation by lifestyle
gg.s.lifestyle <- ggplot(s.df, aes(x=lifestyle, y=S)) +
  geom_violin(colour="grey", lty="dashed", size=0.3) +
  geom_boxplot(aes(fill=lifestyle), width=0.15, size=0.3, outlier.size=0.5) +
  geom_text(data=codonlabelslifestyles.df,
            aes(x=lifestyle, y=Inf, label=test),
            family="mono",
            hjust=0.5,
            vjust=2,
            size=2,
            inherit.aes=FALSE) +
  geom_text(data=plot.labels,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=num),
            fontface="bold",
            hjust=0.5,
            vjust=6,
            size=1.5,
            inherit.aes=FALSE) +
  labs(y="Codon optimisation (S)") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                      labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                    col.df$lifestyle)]),
                      name=NULL) +
  scale_fill_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                    labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                  col.df$lifestyle)]),
                    name=NULL) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(legend.position="none",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.text=element_text(size=6, margin=margin(0, 5, 0, 0)),
        axis.text.x=element_text(colour=c("#009E73","#56B4E9", "#D55E00", "#9AE324", "dimgrey", "#0072B2"),
                                 size=4, face="bold", vjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        panel.grid.major.x=element_blank(),
        plot.margin=margin(10, 10, 15, 10, unit="pt"))

#Format dataframe for plotting pairwise statistical test results
sgenes.pw.df.plot <- sgenes.pw.df %>%
  filter(p.adj.signif != "ns") %>%
  add_xy_position(x="gene_type", scales="free")

#Function to add sample size to strip labels
cust_labeller <- function(x) paste0(x, "\n", plot.labels$num[plot.labels$lifestyle == x])

#Plot boxplots of codon optimisation by gene type by lifestyle
gg.s.genes <- ggplot(s.df2, aes(x=gene_type, y=value)) +
  facet_wrap(~ lifestyle, nrow=1, strip.position="bottom",
             labeller=labeller(lifestyle=as_labeller(cust_labeller))) +
  geom_boxplot(aes(colour=gene_type), position=position_dodge(width=0.5), width=0.4, size=0.3, outlier.size=0.5) +
  stat_pvalue_manual(sgenes.pw.df.plot,
                     label="p.adj.signif", tip.length=0.01, size=1.5, bracket.size=0.2, vjust=-0.05) +
  labs(y="Codon optimisation (S)") +
  scale_colour_grey(start=0.2, end=0.7,
                    labels=c("CSEPs", "CAZymes", "Other"),
                    name=NULL) +
  guides(colour=guide_legend(direction="horizontal")) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(legend.position="top",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.text=element_text(size=6, margin=margin(0, 5, 0, 0)),
        strip.text=element_text(size=4, face="bold", vjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        panel.grid.major.x=element_blank(),
        panel.spacing=unit(0.38, "lines"),
        plot.margin=margin(10, 10, 0, 10, unit="pt"))

#Colour strip labels according to lifestyle
gg.s.genes.grob <- ggplot_gtable(ggplot_build(gg.s.genes))
gg.s.genes.grob.strips <- which(grepl('strip-', gg.s.genes.grob$layout$name))

for (i in seq_along(gg.s.genes.grob.strips)) {
  l <- which(grepl('titleGrob', gg.s.genes.grob$grobs[[gg.s.genes.grob.strips[i]]]$grobs[[1]]$childrenOrder))
  gg.s.genes.grob$grobs[[gg.s.genes.grob.strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <-
    c("#009E73","#56B4E9", "#D55E00", "#9AE324", "dimgrey", "#0072B2")[i]
}

#Write codon optimisation by lifestyle plots to file (Figure 4)
#tiff(file=paste0("Fig4-", Sys.Date(), ".tiff"),
#     height=4, width=3.4, unit="in", res=600, compression="lzw")
plot_grid(gg.s.lifestyle, gg.s.genes.grob, rel_heights=c(0.8, 1), ncol=1, labels="AUTO", label_size=10)
#dev.off()


## LIFESTYLE RANGE SCATTERPLOT ##

#Reformat reported lifestyles
range <- data.frame(lapply(metadata[c((which(colnames(metadata) == "lifestyle") + 1):ncol(metadata))],
                           function(x) str_count(x, paste(sprintf("\\b%s\\b", col.df$lifestyle), collapse = '|'))))

#Make dataframe summarising number of reported lifestyles for each taxon
range.s.df <- data.frame(s.df[c("name", "S")], range=rowSums(range)[match(s.df$taxon, metadata$file2)])
#Remove taxa not identified to species level
range.s.df <- range.s.df[-grep("sp\\.", word(range.s.df$name, 2)),]
#Make second dataframe for grouping strains into one species
range.s.df2 <- data.frame(name=unique(word(range.s.df$name, 2)[which(duplicated(word(range.s.df$name, 2)))]),
                          S=NA,
                          range=NA)

#For each taxon...
for (i in 1:length(range.s.df2$name)) {
  
  #Calculate the mean S value for all grouped strains
  range.s.df2$S[i] <- mean(range.s.df$S[grep(range.s.df2$name[i], range.s.df$name)])
  #Calculate the max number of lifestyles for the species
  range.s.df2$range[i] <- max(range.s.df$range[grep(range.s.df2$name[i], range.s.df$name)])
  
}

#Remove the grouped strains from the original dataframe
range.s.df <- range.s.df[-grep(paste(range.s.df2$name, collapse="|"), range.s.df$name),]
#Combine with the grouped dataframe
range.s.df <- rbind.fill(range.s.df, range.s.df2)

#Make a dataframe for labelling the plot
rangelabels.df <- range.s.df %>% dplyr::count(range)

#Calculate Pearson's correlation
range.corr.r <- signif(cor.test(range.s.df$S, range.s.df$range)$estimate, digits=1)
range.corr.p <- signif(cor.test(range.s.df$S, range.s.df$range)$p.value, digits=1)

#Remove excluded strains from the species tree
plgs.range.tree <- drop.tip(dated.tree.independent$apePhy,
                            s.df$name[union(which(duplicated(word(s.df$name, 2))),
                                            grep("sp\\.", word(s.df$name, 2)))])
#Fix tip labels for the grouped strains
for (i in range.s.df2$name) {
  plgs.range.tree$tip.label[grep(i, plgs.range.tree$tip.label)] <- i
}

#Phylogenetic generalised least squares regression with three different correlation structure methods
range.pgls.pagel <- nlme::gls(S ~ range,
                              correlation=corPagel(1, phy=plgs.range.tree, form=~name),
                              data=range.s.df, method="ML")
range.pgls.brownian <- nlme::gls(S ~ range,
                                 correlation=corBrownian(phy=plgs.range.tree, form=~name),
                                 data=range.s.df, method="ML")
range.pgls.blomberg <- nlme::gls(S ~ range, 
                                 correlation=corBlomberg(1, phy=plgs.range.tree, form=~name, fixed=TRUE),
                                 data=range.s.df, method="ML")

#Select the best model according to AIC 
selected.model <- rownames(anova(range.pgls.pagel, range.pgls.brownian, range.pgls.blomberg))[which.min(anova(range.pgls.pagel, range.pgls.brownian, range.pgls.blomberg)$AIC)]
#Add results for selected model to dataframe
range.s.df$pgls <- predict(get(selected.model))

#Plot scatterplot of range against codon optimisation
gg.range <- ggplot(range.s.df, aes(x=range, y=S)) +
  geom_point(aes(x=as.factor(range)),
             size=0.5, fill="dimgrey", colour="dimgrey") +
  geom_line(aes(y=pgls)) +
  geom_smooth(method="lm", colour="black", linetype="dashed", size=0.5, se=FALSE) +
  geom_text(data=rangelabels.df,
            aes(x=range, y=-Inf, label=paste0("n=", n)),
            fontface="bold",
            colour="dimgrey",
            hjust=0.5,
            vjust=5,
            size=1.5,
            inherit.aes=FALSE) +
  labs(x="Number of reported lifestyles", y="Codon optimisation (S)") +
  annotate("text", x=4, y=0.46, label="PGLS", fontface="bold", size=2) +
  annotate("text", x=4, y=0.44, label=paste0("p=", signif(anova(range.pgls.brownian)[2, 3], 1)), size=2) +
  annotate("text", x=5.5, y=0.46, label="Pearson's\ncorrelation", fontface="bold", size=2) +
  annotate("text", x=5.5, y=0.435, label=paste0("r=", range.corr.r, "\np=", range.corr.p), size=2) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(size=5, face="bold"),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=6, margin=margin(12, 0, 0, 0)),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        plot.margin=margin(10, 10, 0, 10, unit="pt"))

#Write range plot to file (Supplementary Figure 8)
#tiff(file=paste0("SupplementaryFig8-", Sys.Date(), ".tiff"),
#     height=2, width=3.4, unit="in", res=600, compression="lzw")
gg.range
#dev.off()


## CORRELATION WITH PHYLOGENY PCA PLOTS ##

#Read in and format phylogenetic distance matrix from lifestyle comparison test
phyl.dist <- read.csv("lifestyle_comparison/orthogroups/phyldistmatrix.csv")
rownames(phyl.dist) <- phyl.dist$X
phyl.dist$X <- NULL

#Recreate PCA of phylogenetic distances
phylpca <- prcomp(phyl.dist, rank=6)
pca.df <- as.data.frame(phylpca$x)
#Add species complex and S values
pca.df$speciescomplex.abb <- metadata$speciescomplex.abb[match(rownames(pca.df), metadata$short.tip)]
pca.df$s.other <- s.df$S.other[match(metadata$name[match(rownames(pca.df), metadata$short.tip)], s.df$name)]
pca.df$s.CSEP <- s.df$S.CSEP[match(metadata$name[match(rownames(pca.df), metadata$short.tip)], s.df$name)]
pca.df$s.CAZyme <- s.df$S.CAZyme[match(metadata$name[match(rownames(pca.df), metadata$short.tip)], s.df$name)]

#For CSEPs and non-CSEPs...
for (i in c("CSEP", "CAZyme", "other")) {
  
  #Fit S values to PCA
  pca.ordi.s <- ordisurf(as.formula(paste0("phylpca ~ s.", i)), data=pca.df, plot=FALSE)
  
  if (summary(pca.ordi.s)$s.table[4] < 0.05) {
    
    #Pull out coordinates for plotting
    pca.ordi.s.gg <- expand.grid(x=pca.ordi.s$grid$x, y=pca.ordi.s$grid$y)
    #Get z scores
    pca.ordi.s.gg$z <- as.vector(pca.ordi.s$grid$z)
    #Remova NAs
    pca.ordi.s.gg <- data.frame(na.omit(pca.ordi.s.gg))
    
    assign(paste0("pca.ordi.s.", i), pca.ordi.s)
    assign(paste0("pca.ordi.s.gg.", i), pca.ordi.s.gg)
    
  } else {
    
    print(paste0(i, " S values can't be fit to PCA: p=", signif(summary(pca.ordi.s)$s.table[4], 1)))
    
  }
  
}

#Combine S coordinates into one dataframe for plotting
pca.ordi.df <- rbind(data.frame(data="CSEPs",
                                pca.ordi.s.gg.CSEP),
                     data.frame(data="Other",
                                pca.ordi.s.gg.other))

#Calculate centroids for each species complex
pca.centroids.df <- aggregate(. ~ speciescomplex.abb, pca.df[c("PC1", "PC2", "speciescomplex.abb")], mean)

#Make dataframe for labelling plot
pcalabels.df <- data.frame(R=c(signif(summary(pca.ordi.s.CSEP)$r.sq, 1),
                               signif(summary(pca.ordi.s.other)$r.sq, 1)),
                           p=c(signif(summary(pca.ordi.s.CSEP)$s.table[4], 1),
                               signif(summary(pca.ordi.s.other)$s.table[4], 1)),
                           data=c("CSEPs", "Other"))

#Plot PCA of phylogenetic distances with S values fitted
gg.pca <- ggplot(pca.centroids.df,
                 aes(x=PC1, y=PC2,
                     colour=speciescomplex.abb,
                     fill=speciescomplex.abb,
                     shape=speciescomplex.abb)) +
  facet_wrap(~ data, scales="free", ncol=1) +
  geom_contour(data=pca.ordi.df, 
               aes(x=x, y=y, z=z),
               colour="grey",
               binwidth=0.02,
               show.legend=FALSE,
               size=0.3,
               inherit.aes = FALSE) +
  geom_label_contour(data=pca.ordi.df,
                     aes(x=x, y=y, z=z),
                     colour="dimgrey",
                     binwidth=0.04,
                     size=1.5,
                     label.size=NA,
                     label.padding=unit(0.1, "lines"),
                     inherit.aes=FALSE) +
  geom_point(size=1.5, stroke=0.5, position=position_jitter(width=0.05, height=0.05, seed=1)) +
  scale_shape_manual(values=c(1:13)) +
  annotate("text", x=-1.2, y=-1.15, label="ordisurf", fontface="bold", size=1.5) +
  geom_text(data=pcalabels.df,
            aes(x=-1.2, y=-1.5, 
                label=paste0("adj. R=", R,
                             "\np=", p)),
            size=1.5,
            inherit.aes=FALSE) +
  labs(x=paste0("PC1 (", signif(summary(phylpca)$importance[2, 1] * 100, 2), "%)"),
       y=paste0("PC2 (", signif(summary(phylpca)$importance[2, 2] * 100, 2), "%)")) +
  theme_minimal() +
  theme(strip.text=element_text(size=6, face="bold"),
        panel.border=element_rect(colour="black", fill=NA, size=0.5),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_text(size=6),
        aspect.ratio=1,
        panel.grid=element_blank(),
        legend.position="none",
        plot.margin=margin(0, 0, 0, 5))

#Add whether taxa are Fusarium sensu stricto
s.df$genus <- "Allied"
s.df$genus[match(metadata$name[setdiff(grep("Fusarium", metadata$name), grep("'", metadata$name))], s.df$name)] <-
  "Fusarium"

#Test for significant difference in overall S value between Fusarium sensu lato and sensu stricto
#Check for normality of residuals
ggqqplot(residuals(lm(S ~ genus, data=s.df)))
#Check for homogeneity of variances
sgenus.levene <- s.df %>% levene_test(S ~ genus)
sgenus.levene <- data.frame(formula="S ~ genus", sgenus.levene)

#If homogeneity of variance assumption met...
if (sgenus.levene$p > 0.05) {
  
  print("Homogeneity of variance assumption met, doing t-test")
  #Do t-test
  sgenus.ttest <- s.df %>% t_test(S ~ genus)
  #Format
  sgenus.ttest$.y. <- NULL
  sgenus.ttest$formula <- "S ~ genus"
  assign(paste0("sgenus.ttest.", i), sgenus.ttest)
  
} else {
  
  print("Homogeneity of variance assumption not met, doing Wilcoxon rank sum")
  #Do Wilcoxon rank sum
  sgenus.wilcox <- s.df %>% wilcox_test(S ~ genus)
  #Format
  sgenus.wilcox$.y. <- NULL
  sgenus.wilcox$formula <- "S ~ genus"
  assign(paste0("sgenus.wilcox.", i), sgenus.wilcox)
  
}

#Plot boxplot of S values between genera
gg.codongenus <- ggplot(s.df, aes(x=genus, y=S)) +
  geom_violin(colour="grey", lty="dashed", size=0.3) +
  geom_boxplot(width=0.2, size=0.3, outlier.size=0.3) +
  annotate("text", label=paste0("p=", signif(sgenus.ttest$p, 1)), x=1, y=0.57, size=1.2) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text.y=element_text(size=2.5, margin=margin(0, 0, 0, 0)),
        axis.text.x=element_text(size=3.5, face=c("plain", "italic"), margin=margin(0, 0, 0, 0)),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        plot.tag.position=c(1, 0.7),
        plot.margin=margin(3, 3, 3, 3),
        plot.background=element_rect(colour="black", fill="white", size=0.7))

#Remove outgroup from highlighting dataframe
sc.df.pca <- sc.df.spec[!sc.df.spec$sc == "outgroup",]

#Plot tree for PCA legend
gg.pcatree <- ggtree(raxmlng.tree, branch.length="none", linetype=NA) %<+% metadata +
  geom_highlight(data=sc.df.pca, 
                 aes(node=node, fill=sc), extend=-1, alpha=0.5) +
  geom_cladelab(data=sc.df.pca,
                mapping=aes(node=node, label=sc, colour=sc),
                fontsize=1.5, fontface="bold", hjust=0.5, offset.text=6, offset=1.5) +
  geom_tree(size=0.2) +
  scale_x_reverse() +
  coord_cartesian(clip="off") +
  theme(legend.position="none",
        plot.margin=margin(0, 0, 0, 20))

#For each species complex...
for (i in 1:length(sc.df.pca$sc)) {
  
  #Create a dummy plot for the point shape in the PCA
  point <-ggplot(data.frame(x=1,y=1), aes(x, y)) +
    geom_point(shape=which(order(sc.df.pca$sc) == i), size=0.9, stroke=0.3,
               col=hue_pal()(13)[which(order(sc.df.pca$sc) == i)]) +
    scale_x_continuous(breaks=1) +
    theme_void() +
    theme(aspect.ratio=1)
  #Convert to grob
  point.grob <- as.grob(point)
  #Calculate position in tree
  pos <- mean(gg.pcatree$data$y[which(gg.pcatree$data$speciescomplex.abb == sc.df.pca$sc[i])])
  
  #Add point to tree plot
  gg.pcatree <- gg.pcatree +
    annotation_custom(point.grob, xmin=-31, xmax=-33,
                      ymin=pos-1,
                      ymax=pos+1)
}

#Make dataframe for grid of PCs 1-6
pca.grid.df <- pca.df
colnames(pca.grid.df)[1:6] <- paste0("PC", seq(1:6), "\n(",
                                     signif(summary(phylpca)$importance[2, 1:6] * 100, 2), "%)")

#Plot grid of PCs 1-6
gg.pca.grid <- 
  ggplot(pca.grid.df, aes(x=.panel_x, y=.panel_y)) + 
  geom_point(aes(colour=speciescomplex.abb, fill=speciescomplex.abb, shape=speciescomplex.abb), 
             size=0.9, stroke=0.3) + 
  facet_matrix(vars(1:6),
               switch="both",
               flip.rows=TRUE,
               layer.diag=FALSE,
               layer.upper=FALSE) +
  scale_shape_manual(values=c(1:13)) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(strip.text.y=element_text(size=4, face="bold"),
        strip.text.x=element_text(size=4, face="bold", margin=margin(20, 0, 0, 0)),
        strip.placement="outside",
        panel.border=element_rect(colour="black", fill=NA, size=0.4),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_text(size=6),
        aspect.ratio=1,
        panel.grid=element_blank(),
        legend.position="none",
        plot.margin=margin(5, -50, -55, 0))

#Remove empty grid positions
gg.pca.grid.grob <- ggplotGrob(gg.pca.grid)
pos <- gg.pca.grid.grob$layout$name %in% c(paste0("panel-6-", c(1:6)),
                                           paste0("panel-5-", c(2:6)),
                                           paste0("panel-4-", c(3:6)),
                                           paste0("panel-3-", c(4:6)),
                                           paste0("panel-2-", c(5:6)),
                                           paste0("panel-1-", 6), "strip-l-6", "strip-b-6")
gg.pca.grid.grob$grobs <- gg.pca.grid.grob$grobs[!pos]
gg.pca.grid.grob$layout <- gg.pca.grid.grob$layout[!pos, ]

gg.pca.grid.grob$layout$t[gg.pca.grid.grob$layout$name %in% paste0("strip-b-", c(1:5))] <- 16
gg.pca.grid.grob$layout$b[gg.pca.grid.grob$layout$name %in% paste0("strip-b-", c(1:5))] <- 16

#Write grid of PCAs to file (Supplementary Fig 13)
#tiff(file=paste0("SupplementaryFig13-", Sys.Date(), ".tiff"),
#     height=4.8, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(gg.pca.grid.grob, gg.pcatree,
          align="horizontal", axis="t", rel_widths=c(1, 0.4))
#dev.off()


## CODON USAGE BIAS HIERARCHICAL CLUSTERING ##

#Combine RSCU results for all taxa
rscu.grid <- t(cbind(as.data.frame(mget(ls(pattern="rscu.")))))
rownames(rscu.grid) <- metadata$name[match(sub("rscu.", "", rownames(rscu.grid)), metadata$file)]
#Remove outgroup
rscu.grid <- rscu.grid[match(metadata$name[metadata$ingroup == 1], rownames(rscu.grid)),]
#Remove non redundant codons
rscu.grid <- rscu.grid[, !(colnames(rscu.grid) %in% c("tga", "tgg", "tag", "taa", "atg"))]

#Normalise data
rscu.grid <- scale(rscu.grid)
#Make distance matrix
rscu.dist <- dist(rscu.grid, method="euclidean")
#Do hierarchical clustering
rscu.hclust <- hclust(rscu.dist, method="average")

#Calculate topological distance from species tree
RF.dist(as.phylo(rscu.hclust, directed=FALSE), drop.tip(raxmlng.tree, outgroup), normalize=TRUE)
#Number of bipartitions that differ
print(paste0(RF.dist(as.phylo(rscu.hclust, directed=FALSE), drop.tip(raxmlng.tree, outgroup), normalize=TRUE) *
  length(bitsplits(unroot(drop.tip(raxmlng.tree, outgroup)))[3]$freq), "/",
       length(bitsplits(unroot(drop.tip(raxmlng.tree, outgroup)))[3]$freq)))

#Create random trees to compare topological distance to
random.trees <- rep(NA, 1000)

for (i in 1:length(random.trees)){
  
  r.tree <- rtree(length(labels(rscu.hclust)), rooted=FALSE, tip.label=labels(rscu.hclust), br=NULL)
  
  random.trees[i] <- RF.dist(r.tree, drop.tip(raxmlng.tree, outgroup), normalize=TRUE)
  
}

#p value for hierarchical clustering being similar to species tree
sum(random.trees < RF.dist(as.phylo(rscu.hclust, directed=FALSE),
                           drop.tip(raxmlng.tree, outgroup), normalize=TRUE)) /
  length(random.trees)
  
#Plot hierarchical clustering
gg.hclust <- ggtree(rscu.hclust, size=0.2) %<+% metadata +
  xlim(0,20) +
  geom_tiplab(aes(label=tiplab2), offset=7, size=1, face="italic") +
  geom_tippoint(aes(col=lifestyle), size=0.7) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                      labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                    col.df$lifestyle)]),
                      na.translate=FALSE,
                      guide=guide_legend(title="Lifestyle",
                                         ncol=2,
                                         title.hjust=0,
                                         order=1)) +
  theme(legend.box="vertical",
        legend.title=element_text(size=5, face="bold"),
        legend.text=element_text(size=5, margin=margin(0, 3, 0, 0)),
        legend.key.size=unit(6, "pt"),
        legend.spacing=unit(0, "pt"),
        plot.margin=margin(0, 0, 0, 0, unit="pt"))

#Add heatmap of RSCU values
gg.rscu <- gheatmap(gg.hclust,
                    width=0.7,
                    as.data.frame(rscu.grid),
                    color="grey",
                    offset=0.2,
                    colnames=FALSE) +
  scale_fill_gradient2(low="#F0E442", mid="white", high="#CC79A7",
                       breaks=pretty_breaks(),
                       guide=guide_colourbar(title="Normalised\ncodon usage bias\n(RSCU)",
                                             title.position="top",
                                             direction="horizontal",)) +
  theme(legend.position=c(0.2, 0.8))

## COMBINE ##

#Write to file
#tiff(file=paste0("Fig5-", Sys.Date(), ".tiff"),
#     height=4, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.pca,
                    gg.pcatree +
                      annotation_custom(ggplotGrob(gg.codongenus), ymin=55, ymax=65, xmin=-25, xmax=-45),
                    rel_widths=c(1.5, 1), align="h", axis="tb",
                    labels=c("A", ""), label_size=10),
          gg.rscu,
          labels="AUTO", label_size=10, ncol=2, rel_widths=c(0.85, 1))
#dev.off()


#####################    SUPPLEMENTARY FIGURE 9    #############################
#################   CSEP PREDICTION PIPELINE SANKEY   ##########################

#Make empty list to capture protein lengths
protein.lengths <- list()

#For each taxon...
for (i in metadata$file2) {
  
  #Read in and format CSEPfilter log
  CSEPfilter <- read.csv(paste0("CSEP_CAZyme_prediction/CSEPfilter_", i, ".faa.log"),
                         sep=":", row.names=NULL, header=FALSE)
  CSEPfilter <- CSEPfilter[!is.na(CSEPfilter$V2),]
  CSEPfilter$V1 <- word(CSEPfilter$V1, 1)
  CSEPfilter$V1[max(grep("Phobius", CSEPfilter$V1))] <- "Phobius2"
  rownames(CSEPfilter) <- CSEPfilter$V1
  CSEPfilter <- subset(CSEPfilter, select="V2")
  colnames(CSEPfilter) <- metadata$name[metadata$file2 == i]
  
  assign(paste0(i, ".CSEPfilter"), CSEPfilter)
  
  #Read in CSEPs
  CSEPs <- scan(paste0("CSEP_CAZyme_prediction/", i, ".faa_candidate_effectors"), character(), quote="")
  
  assign(paste0(i, ".CSEPs"), CSEPs)
  
  #Read in proteins
  proteins <- read.fasta(paste0("orthology_inference/", i, ".faa"), seqtype="AA")
  names(proteins) <- sub(".*\\.faa_", "", names(proteins))
  
  #Create bar to show progress
  progress.bar <- txtProgressBar(1, length(proteins), initial=0, char="=", style=3)
  
  counter <- 0
  
  #For each protein...
  for (j in names(proteins)) {
    
    counter <- counter + 1
    
    #Update progress bar
    setTxtProgressBar(progress.bar, counter)
    
    #Get length
    protein.lengths[[i]][j] <- length(proteins[[j]])
    
  }
  
}

#Combine CSEPfilter results into master dataframe
CSEPfilter.df <- bind_cols(mget(paste0(metadata$file2, ".CSEPfilter")))

#Make dataframe for labelling plot
CSEPfilterlabels.df <- data.frame(programme=c("", "SignalP v5.0b", "TargetP v2.0", "Phobius v1.01", "TMHMM v2.0c", "Phobius v1.01", "ps_scan v1.86", "NucPred v1.1",  "PredGPI", "EffectorP v3.0"),
                                  x=rev(seq(1.55, length(rownames(CSEPfilter.df)) + 1)),
                                  mean=round(rowMeans(CSEPfilter.df)),
                                  range=paste0(prettyNum(rowMins(as.matrix(CSEPfilter.df)), big.mark=","),
                                               " - ", 
                                               prettyNum(rowMaxs(as.matrix(CSEPfilter.df)), big.mark=",")),
                                  sum=rowSums(CSEPfilter.df[,-grep("programme", colnames(CSEPfilter.df))]),
                                  arrow.x=c(seq(2, 10), NA),
                                  arrow.xend=c(seq(1.2, 9.2), NA),
                                  face=c("bold", "plain", "plain", "bold", "plain",
                                         "plain", "plain", "plain", "plain", "bold"))

#Add column with programme
CSEPfilter.df$programme <- rownames(CSEPfilter.df)
#Melt dataframe for plotting
CSEPfilter.df <- melt(CSEPfilter.df)
#Order programmes
CSEPfilter.df$programme <- factor(CSEPfilter.df$programme,
                                  levels=c("Total", "SignalP", "TargetP", "Phobius", "TMHMM",
                                           "Phobius2", "Prosite", "NucPred", "PredGPI", "EffectorP"))

#Plot sankey diagram of CSEPfilter steps
gg.CSEPfilter <- ggplot(CSEPfilter.df,
                        aes(x=programme, y=value, stratum=variable, alluvium=variable, fill=variable)) +
  geom_stratum(colour=NA, size=0.2) +
  geom_flow(size=0.2) +
  geom_segment(data=CSEPfilterlabels.df,
               aes(x=arrow.x, xend=arrow.xend, y=-130000, yend=-130000),
               arrow=arrow(length=unit(0.1, "cm"), type="closed"),
               size=0.15,
               inherit.aes=FALSE) +
  geom_label(data=CSEPfilterlabels.df,
             aes(x=rev(c(1:length(rownames(CSEPfilterlabels.df)))), y=-130000, label=range),
             size=2,
             inherit.aes=FALSE) +
  geom_label(data=CSEPfilterlabels.df[CSEPfilterlabels.df$programme != "",],
             aes(x=x, y=-130000, label=programme),
             label.padding=unit(0.1, "lines"),
             label.size=NA,
             size=1.5,
             inherit.aes=FALSE) +
  annotate("text", x=10.5, y=-130000, label="Number of proteins", size=2) +
  labs(tag="A") +
  scale_x_discrete(limits=rev(c("Blank",
                                "Total",
                                "SignalP",
                                "TargetP",
                                "Phobius",
                                "TMHMM",
                                "Phobius2",
                                "Prosite",
                                "NucPred",
                                "PredGPI",
                                "EffectorP")),
                   labels=rev(c("",
                                "Predicted\nproteomes",
                                "Predicted\nsignal peptide",
                                "",
                                "Predicted\nsecretomes",
                                "TM helices\n<=1",
                                "",
                                "No ER retention",
                                "Not nuclear\nlocalised",
                                "Not\nGPI-anchored",
                                "CSEPs")),
                   expand=c(0, 0)) +
  scale_fill_manual(values=rep(c("dimgrey", "grey"), length.out=length(metadata$file2))) +
  coord_flip(ylim=c(-200000, sum(CSEPfilter.df$value[CSEPfilter.df$programme == "Total"])),
             clip="off") +
  theme(legend.position="none",
        plot.tag=element_text(face="bold", size=10),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=8, face=rev(CSEPfilterlabels.df$face)),
        axis.text.x=element_blank(),
        axis.title=element_blank())

#Function to add brackets
#https://stackoverflow.com/questions/35633239/add-curly-braces-to-ggplot2-and-then-use-ggsave
bracketsGrob <- function(...){
  l <- list(...)
  e <- new.env()
  e$l <- l
  grid:::recordGrob(  {
    do.call(grid.brackets, l)
  }, e)
}

#Plot brackets
b1 <- bracketsGrob(0.025, 0.68, 0.025, 0.83, lwd=1.5, ticks=0.8, col="dimgrey", h=0.02)
b2 <- bracketsGrob(0.025, 0.38, 0.025, 0.53, lwd=1.5, ticks=0.8, col="dimgrey", h=0.02)

#Add brackets to sankey
gg.CSEPfilter <- gg.CSEPfilter +
  annotation_custom(b1) +
  annotation_custom(b2)

#Make dataframe of protein lengths
lengths.df <- stack(protein.lengths)
#Make vector with CSEPs from all taxa
CSEPs <- unlist(mget(paste0(metadata$file2, ".CSEPs")))
#Label proteins as CSEPs or not
lengths.df$group[!is.na(match(rownames(lengths.df), CSEPs))] <- "CSEP"
lengths.df$group[is.na(lengths.df$group)] <- "other"

#Make dataframe for labelling plot
lengthslabels.df <- data.frame(table(lengths.df$group))
lengthslabels.df$mean <- t.test(lengths.df$values[lengths.df$group == "CSEP"], lengths.df$values[lengths.df$group == "other"])$estimate
lengthslabels.df$label <- paste0("n=", prettyNum(lengthslabels.df$Freq, big.mark=","))

#Plot boxplot of protein lengths
gg.lengths <- ggplot(lengths.df, aes(x=group, y=values)) +
  geom_hline(yintercept=300, lty="dashed") +
  geom_violin(size=0.3) +
  geom_boxplot(width=0.2, size=0.3, outlier.size=0.3) +
  geom_text(data=lengthslabels.df, aes(x=Var1, y=mean, label=label),
            nudge_y=500, nudge_x=0.3, size=1.5, inherit.aes=FALSE) +
  scale_x_discrete(labels=c("CSEP", "Other"),
                   limits=c("CSEP", "other")) +
  scale_y_continuous(labels=addUnits, 
                     sec.axis=dup_axis(),
                     breaks=c(300, 5000, 10000),
                     minor_breaks=seq(0, round_any(max(lengths.df$values), 1000, f=ceiling), 1000)) +
  labs(y="Protein length (aa)", tag="B") +
  theme_minimal() +
  theme(panel.grid.major.x=element_blank(),
        plot.tag=element_text(face="bold", size=10),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.y.left=element_text(size=10),
        axis.title.y.right=element_blank(),
        axis.text.y.left=element_blank(),
        axis.text.y.right=element_text(size=8),
        panel.border=element_rect(colour="black", fill=NA, size=1.5))

#Write to file
#tiff(file=paste0("SupplementaryFig9-", Sys.Date(), ".tiff"),
#     height=4, width=6.75, units="in", res=600, compression="lzw")
gg.CSEPfilter + 
  annotation_custom(ggplotGrob(gg.lengths), ymin=300000, ymax=900000, xmin=2, xmax=8)
#dev.off()



## SIMPLIFIED GRAPHICAL ABSTRACT TREE

#Make dataframe for plotting grid of lifestyles on circular tree
lifestyles.df <- metadata[,c(which(colnames(metadata) == "name"),
                             c((which(colnames(metadata) == "lifestyle")+1):length(colnames(metadata))))]
lifestyles.df <- lifestyles.df %>% mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
lifestyles.df[2:length(colnames(lifestyles.df))][!is.na(lifestyles.df[2:length(colnames(lifestyles.df))])] <- 1
lifestyles.df[2:length(colnames(lifestyles.df))][is.na(lifestyles.df[2:length(colnames(lifestyles.df))])] <- 0
lifestyles.df <- melt(lifestyles.df, id="name")

library(ggtreeExtra)

#Plot tree
gg.abstractree <- ggtree(dated.tree.independent$apePhy, linetype=NA, layout="circular") %<+% metadata +
  geom_highlight(data=sc.df.dated, 
                 aes(node=node, fill=as.factor(box)), alpha=1, extend=2, show.legend=FALSE) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  new_scale_fill() +
  geom_tree(size=0.25) +
  geom_fruit(data=lifestyles.df, geom=geom_tile,
             mapping=aes(y=name, x=variable, alpha=value, fill=variable),
             color="grey50", offset=0.05,size=0.02) +
  scale_alpha_discrete(range=c(0, 1), guide=FALSE) +
  geom_tippoint(aes(colour=lifestyle), size=1, show.legend=FALSE) +
  scale_fill_manual(values=col.df$colour[match(levels(col.df$lifestyle), col.df$lifestyle)],
                    labels=str_to_sentence(col.df$lifestyle[match(levels(col.df$lifestyle), col.df$lifestyle)])) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(metadata$lifestyle[match(dated.tree$apePhy$tip.label,
                                                                                      metadata$name)])),
                                                 col.df$lifestyle)]) +
  guides(fill=guide_legend(direction="horizontal",
                           title=NULL,
                           ncol=2)) +
  theme(legend.position=c(0.51, 0.43),
        legend.background=element_rect(fill=NA),
        legend.key.size=unit(0.1, "cm"),
        legend.text=element_text(size=5, margin=margin(0,0,2,-2)),
        plot.margin=margin(0, 0, 0, 0))

gg.abstractree

#tiff(file=paste0("GraphicAbstractTree-", Sys.Date(), ".tiff"),
#     height=4, width=4, unit="in", res=300, compression="lzw")
gg.abstractree
#dev.off()


## MISC STATS FOR TEXT

#Num of both CSEPs and CAZymes
length(intersect(which(!is.na(orthogroups.stats.ingroup0$CSEP)), which(!is.na(orthogroups.stats.ingroup0$CAZyme))))

#Copynum range of PL1_7
range(cazyme.fam.df$value[cazyme.fam.df$variable == "PL1_7"])

#Mutant phenotypes of positively selected genes
phibase.df[match(unique(na.omit(absrel.df$CSEP_name)), phibase.df$Gene), c("Gene", "Gene.Function", "Pathogen.species", "Host.species", "Mutant.Phenotype")]

#Names of positively selected CAZymes
orthogroups.stats.ingroup0[match(paste0("OG000", unique(unlist(strsplit(absrel.df$CAZyme[which(!is.na(absrel.df$CAZyme))], " ")))), orthogroups.stats.ingroup0$orthogroup), c("orthogroup", "CAZyme_family", "CAZyme", "EC")]
#Known substrates for positively selected CAZYmes
substrates.df[substrates.df$CAZy.Family %in% orthogroups.stats.ingroup0$CAZyme_family[match(paste0("OG000", unique(unlist(strsplit(absrel.df$CAZyme[which(!is.na(absrel.df$CAZyme))], " ")))), orthogroups.stats.ingroup0$orthogroup)],]


## Write tables summarising all statistical test results

levene.table <- rbind(slifestyles.levene,
                      do.call("rbind", mget(ls(pattern="\\.levene"))))

levene.table$p <- signif(levene.table$p, 1)
levene.table$formula <- sub("\n", " ", levene.table$formula)

anova.table <- as.data.frame(rbind(num.genes.anova.orthogroup,
                                   num.genes.anova.CSEP,
                                   num.genes.anova.CAZyme,
                                   
                                   specific.anova.Orthogroups,
                                   specific.anova.CSEPs,
                                   
                                   copynum.anova.Orthogroups,
                                   copynum.anova.CSEPs,
                                   copynum.anova.CAZymes,
                                   
                                   contrastfel.anova.increase$significance,
                                   contrastfel.anova.decrease,
                                   
                                   slifestyles.anova,
                                   
                                   sgenes.anova.endophyte$significance,
                                   sgenes.anova.insectmutualist,
                                   sgenes.anova.plantassociate,
                                   sgenes.anova.plantpathogen$significance,
                                   sgenes.anova.saprotroph))

anova.table$p <- signif(anova.table$p, 1)
anova.table$formula <- gsub("\n", " ", anova.table$formula)

pw.table <- rbind(slifestyles.pw.df[-c(1, 4, 9)], 
                  sgenes.pw.df[-c(1, 2, 9)],
                  contrastfel.pw.df.increase)

pw.table$p.adj <- signif(pw.table$p.adj, 1)
pw.table$formula <- gsub("\n", " ", pw.table$formula)

permanova.table <- 
  rbind(data.frame(group="Orthogroups", effect=rownames(permanova.orthogroups), permanova.orthogroups),
        data.frame(group="CSEPs", effect=rownames(permanova.CSEPs), permanova.CSEPs),
        data.frame(group="CAZymes", effect=rownames(permanova.CAZymes), permanova.CAZymes))

#write.csv(levene.table, "levene_results.csv", row.names=FALSE)
#write.csv(anova.table, "anova_results.csv", row.names=FALSE)
#write.csv(pw.table, "pairwise_results.csv", row.names=FALSE)
#write.csv(permanova.table, "permanova_results.csv", row.names=FALSE)
