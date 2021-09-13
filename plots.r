###################################
###################################
####                           ####
####  Script to plot figures   ####
####                           ####
###################################
###################################

library(ape)
library(cowplot)
library(deeptime)
library(dendextend)
library(plyr)
library(dplyr)
library(eulerr)
library(ggplot2)
library(ggalluvial)
library(ggnewscale)
library(ggplotify)
library(ggpubr)
library(ggrepel)
library(ggthemes)
library(ggtree)
library(grid)
library(jsonlite)
library(matrixStats)
library(MCMCtreeR)
library(metR)
library(multcompView)
library(nlme)
library(pBrackets)
library(phytools)
library(reshape2)
library(scales)
library(seqinr)
library(stringr)
library(vegan)

#Colour palette
show_col(colorblind_pal()(8))

#Read in orthogroup data
load("CSEP_prediction/orthogroup-matrices-2021-09-09.RData")

#Read in sample metadata
metadata <- read.csv("metadata.csv")

#Make dataframe of lifestyle colours
col.df <- data.frame(lifestyle=c("endophyte", "animal pathogen", "human pathogen","animal associate",
                                 "insect mutualist", "plant associate", "plant pathogen", "saprotroph",
                                 "mycoparasite"),
                     colour=c("#009E73", "#FFE983", "#000000", "#F1BCF4", "#56B4E9",
                              "#9AE324", "dimgrey", "#0072B2", "#D55E00"))


############################    FIGURE 2    ####################################
#################   CSEP PREDICTION PIPELINE SANKEY   ##########################

#Make empty list to capture protein lengths
protein.lengths <- list()

#For each taxon...
for (i in metadata$file2) {
  
  #Read in and format CSEPfilter log
  CSEPfilter <- read.csv(paste0("CSEP_prediction/CSEPfilter_", i, ".faa.log"),
                         sep=":", row.names=NULL, header=FALSE)
  CSEPfilter <- CSEPfilter[!is.na(CSEPfilter$V2),]
  CSEPfilter$V1 <- word(CSEPfilter$V1, 1)
  CSEPfilter$V1[max(grep("Phobius", CSEPfilter$V1))] <- "Phobius2"
  rownames(CSEPfilter) <-CSEPfilter$V1
  CSEPfilter <- subset(CSEPfilter, select="V2")
  colnames(CSEPfilter) <- metadata$name[metadata$file2 == i]
  
  assign(paste0(i, ".CSEPfilter"), CSEPfilter)
  
  #Read in CSEPs
  CSEPs <- scan(paste0("CSEP_prediction/", i, ".faa_candidate_effectors"), character(), quote="")
  
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
        plot.tag=element_text(face="bold"),
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
  scale_y_continuous(labels=comma, sec.axis=dup_axis(),
                     breaks=c(300, 5000, 10000),
                     #expand=c(0, 100),
                     minor_breaks=seq(0, round_any(max(lengths.df$values), 1000, f=ceiling), 1000)) +
  labs(y="Protein length (aa)", tag="B") +
  theme_minimal() +
  theme(panel.grid.major.x=element_blank(),
        plot.tag=element_text(face="bold"),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.y.left=element_text(size=10),
        axis.title.y.right=element_blank(),
        axis.text.y.left=element_blank(),
        axis.text.y.right=element_text(size=8),
        panel.border=element_rect(colour="black", fill=NA, size=1.5))

#Write to file
tiff(file=paste0("Fig2-", Sys.Date(), ".tiff"),
     height=4, width=6.75, units="in", res=600, compression="lzw")
gg.CSEPfilter + 
  annotation_custom(ggplotGrob(gg.lengths), ymin=300000, ymax=900000, xmin=2, xmax=8)
dev.off()


############################    FIGURE 3    ####################################
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
                        sep="\t", skip=644, nrows=61, header=FALSE)
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

#Write to file (Supplementary Figure 1)
tiff(file=paste0("SupplementaryFig1-", Sys.Date(), ".tiff"),
     height=4, width=6.75, units="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.mcmc.correlated, gg.infinitesite.correlated, align="h", axis="tb"),
          plot_grid(gg.mcmc.independent, gg.infinitesite.independent, align="h", axis="tb"), nrow=2,labels="AUTO")
dev.off()


## DATED PHYLOGENY ##

#Read in trees
astral <- read.tree("phylogenomics/species_tree/astral/fus_astral_proteins_62T.tre")
astral$edge.length <- rep(1, length(astral$edge.length))
raxmlng <- read.tree("phylogenomics/species_tree/raxml-ng/fus_proteins_62T.raxml.support")
iqtree <- read.tree("phylogenomics/species_tree/iqtree/fus_proteins_62T_iqtree_genepart.contree")
#Make vector with outgroup
outgroup <- "Ilyonectria sp."

#For each species-tree method...
for (i in c("iqtree", "raxmlng", "astral")) {
  
  #Get the tree
  tree <- get(i)
  #Edit tip names
  tree$tip.label <- metadata$name[match(tree$tip.label, metadata$tip)]
  #Root tree
  tree <- root(tree, outgroup, resolve.root=TRUE, edgelabel=TRUE)
  #Ultrametrise
  chrono <- suppressWarnings(chronos(tree))
  #Convert to dendrogram for tanglegrams
  dend <- as.dendrogram(chrono)
  
  assign(paste0(i, ".dend"), dend)
  
  #Mark unsupported branches
  if (i == "iqtree") {
    tree$node.label[which(suppressWarnings(as.numeric(tree$node.label) < 95))] <- "∗"
    tree$node.label[tree$node.label != "∗"] <- ""
  } else if (i == "raxmlng") {
    tree$node.label[which(suppressWarnings(as.numeric(tree$node.label) < 70))] <- "†"
    tree$node.label[tree$node.label != "†"] <- ""
  } else if (i == "astral") {
    tree$node.label[which(suppressWarnings(as.numeric(tree$node.label) < 0.95))] <- "×"
    tree$node.label[tree$node.label != "×"] <- ""
  }
  
  #Plot tree
  gg <- ggtree(tree, branch.length="none") %<+% metadata
  #Capture data structure
  gg.tree.data <- gg[["data"]] %>%
    arrange(y)
  
  assign(paste0(i, ".tree"), tree)
  assign(paste0("gg.", i), gg)
  assign(paste0("gg.tree.data.", i), gg.tree.data)
  
}

#Untangle dendrograms and plot to compare topologies
dend <- untangle_labels(raxmlng.dend, iqtree.dend, method="step2side")
tanglegram(dend,
           main_left="Concatenated ML (RAxML-NG)",
           main_right="Coalescent (ASTRAL-III)",
           lab.cex=1,
           cex_main=1,
           columns_width=c(2,0.5,2),
           lwd=1,
           margin_outer=0,
           margin_inner=16,
           margin_top=1.5,
           margin_bottom=0,
           axes=FALSE,
           edge.lwd=1,
           color_lines="grey",
           rank_branches=TRUE)

dend <- untangle_labels(raxmlng.dend, astral.dend, method="step2side")
tanglegram(dend,
           main_left="Concatenated ML (RAxML-NG)",
           main_right="Coalescent (ASTRAL-III)",
           lab.cex=1,
           cex_main=1,
           columns_width=c(2,0.5,2),
           lwd=1,
           margin_outer=0,
           margin_inner=16,
           margin_top=1.5,
           margin_bottom=0,
           axes=FALSE,
           edge.lwd=1,
           color_lines="grey",
           rank_branches=TRUE)

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

#Get matching node numbers between trees
ml.nodes <- matchNodes(raxmlng.tree, iqtree.tree)
ml.nodes2 <- matchNodes(astral.tree, iqtree.tree)

#Assign species tree topology
species.tree <- iqtree.tree
#Combine branch supports from all tree building methods on the species tree
species.tree$node.label <- paste0(species.tree$node.label,
                                  gg.tree.data.raxmlng$label[match(ml.nodes[,1],
                                                                   gg.tree.data.raxmlng$node)][order(ml.nodes[,2])])
species.tree$node.label <- paste0(species.tree$node.label,
                                  gg.tree.data.astral$label[match(ml.nodes2[,1],
                                                                  gg.tree.data.astral$node)][order(ml.nodes2[,2])])

#Plot species tree
gg.speciestree <- ggtree(species.tree, branch.length="none") %<+% metadata

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
  nodes <- matchNodes(species.tree, dated.tree$apePhy)
  dated.tree$apePhy$node.label <- gg.tree.data.speciestree$label[
    match(nodes[,1], gg.tree.data.speciestree$node)][order(nodes[,2])]
  
  assign(paste0("dated.tree.", clock), dated.tree)
  
  #Plot tree
  gg.datedtree <-  ggtree(dated.tree$apePhy, linetype=NA) %<+% metadata
  
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
    theme(plot.margin=unit(c(1, 0, 0, 0), "cm"),
          axis.text.x.bottom=element_text(size=5),
          axis.title.x.bottom=element_text(size=5))
  
  #Reverse the x axis
  gg.datedtree <- revts(gg.datedtree)

  #Make dataframe of species complex nodes
  sc.df.dated <- data.frame(sc=names(table(metadata$speciescomplex.abb[
    match(dated.tree$apePhy$tip.label, metadata$name)]))[
      table(metadata$speciescomplex.abb[match(dated.tree$apePhy$tip.label, metadata$name)]) > 1],
                            node=NA)
  
  #Get nodes for each species complex
  for (i in 1:length(sc.df.dated$sc)) {
    sc.df.dated$node[i] <- getMRCA(dated.tree$apePhy,
                                   metadata$name[metadata$speciescomplex.abb == sc.df.dated$sc[i]])
  }
  
  #Make dataframe for species complexes with only one taxon represented
  sc.df.dated2 <- data.frame(sc=unique(metadata$speciescomplex.abb)[
    is.na(match(unique(metadata$speciescomplex.abb), sc.df.dated$sc))],
                             node=NA)
  
  sc.df.dated2$node <- gg.tree.data.dated$node[match(sc.df.dated2$sc, gg.tree.data.dated$speciescomplex.abb)]
  
  #Combine dataframes
  sc.df.dated <- rbind(sc.df.dated, sc.df.dated2)
  sc.df.dated <- sc.df.dated[match(na.omit(unique(gg.tree.data.dated$speciescomplex.abb)), sc.df.dated$sc),]
  #Make alternated coding for highlights on tree
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
    geom_nodelab(size=1.4, nudge_x=-0.5, hjust=1, vjust=-0.5) +
    geom_tiplab(aes(label=tiplab), fontface=tiplabel.face.datedtree, size=1.5, offset=1) +
    geom_tree(size=0.25) +
    geom_tippoint(aes(colour=lifestyle), size=0.7, show.legend=FALSE) +
    scale_colour_manual(values=col.df$colour[match(sort(unique(metadata$lifestyle[match(dated.tree$apePhy$tip.label,
                                                                                        metadata$name)])),
                                                   col.df$lifestyle)]) +
    scale_linetype_manual(values=c("solid", "11"))
  
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

#Write dated trees with confidence interval tables to file (Supplementary Figure 4)
tiff(file=paste0("SupplementaryFig4-", Sys.Date(), ".tiff"),
     height=8, width=6.75, unit="in", res=600, compression="lzw")
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
                      theme(plot.margin=margin(0, 1, 0, 0, unit="cm")),
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
                      theme(plot.margin=margin(0, 1, 0, 0, unit="cm")),
                    gg.ci.first.independent, gg.ci.second.independent,
                    rel_widths=c(6, 1 ,1), nrow=1, align="h", axis="tb"),
          ncol=1, labels="AUTO", align="v", axis="lr")
dev.off()


## NUMBER OF GENES/CSEPs BARGRAPH ##

#Make dataframe categorising genes into core, accessory or specific
genes.df <- data.frame(taxon=rep(colnames(orthogroups.copies.ingroup1), each=3), 
                       category=rep(c("specific", "accessory", "core"),
                                    length(colnames(orthogroups.copies.ingroup1))),
                       orthogroups=NA,
                       CSEPs=NA)

#Filter orthogroup presence-absence matrix for orthogroups with at least one predicted CSEP
CSEP.orthogroups <- orthogroups.copies.ingroup1[
  match(rownames(CSEP.count.ingroup1)[which(rowSums(CSEP.count.ingroup1) > 0)],
        rownames(orthogroups.copies.ingroup1)),]

#Check for correlation between N50 and number of CSEPs
CSEPcheck.df <- data.frame(CSEPs=colSums(CSEP.orthogroups))
CSEPcheck.df$taxon <- rownames(CSEPcheck.df)
CSEPcheck.df$N50 <- metadata$N50[match(CSEPcheck.df$taxon, metadata$file2)]
cor.test(CSEPcheck.df$CSEPs, CSEPcheck.df$N50)

#For each taxon...
for (i in unique(genes.df$taxon)) {
  
  #For each category...
  for (j in unique(genes.df$category)) {
    
    #Get number of orthogroups
    genes.df$orthogroups[intersect(which(genes.df$taxon == i), which(genes.df$category == j))] <-
      table(orthogroups.stats.ingroup1$category[
        match(rownames(orthogroups.copies.ingroup1[orthogroups.copies.ingroup1[, i] > 0,]),
              orthogroups.stats.ingroup1$orthogroup)])[j]
    
    #Get number of CSEPs
    genes.df$CSEPs[intersect(which(genes.df$taxon == i), which(genes.df$category == j))] <-
      table(orthogroups.stats.ingroup1$category[
        match(rownames(CSEP.orthogroups[CSEP.orthogroups[, i] > 0,]),
              orthogroups.stats.ingroup1$orthogroup)])[j]
    
  }
}

#Add taxon name
genes.df$name <- metadata$name[match(genes.df$taxon, metadata$file2)]
#Order categories
genes.df$category <- factor(genes.df$category, levels=c("specific", "accessory", "core"))
#Order taxa by tips in the tree
genes.df$name <- factor(genes.df$name, levels=rev(tip.order.datedtree))
#Add column for order of taxa in tree
genes.df$y <- as.numeric(factor(genes.df$name))

#Plot bargraph for number of orthogroups
gg.category.orthogroups <- ggplot(genes.df, aes(y=name, x=orthogroups, fill=category)) +
  geom_bar(stat="identity", size=0.5, width=0.6) +
  labs(x="All genes") +
  scale_y_discrete(limits=rev(tip.order.datedtree),
                   expand=expansion(0, 0)) +
  coord_cartesian(clip="off",
                  ylim=c(0, 62)) +
  guides(fill=guide_legend(nrow=1)) +
  scale_x_continuous(expand=c(0, 0),
                     position="top",
                     labels=scales::comma) +
  scale_fill_manual(values=c("lightgrey", "darkgrey", "dimgrey"),
                    breaks=c("core", "accessory", "specific"),
                    labels=c("Core", "Accessory", "Specific")) +
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x.top=element_text(size=6, margin=margin(0, 0, 4, 0)),
        axis.text.x.top=element_text(size=3),
        panel.grid.major.x=element_line(colour="white"),
        panel.grid.minor.x=element_line(colour="white"),
        panel.grid.major.y=element_blank(),
        legend.position=c(0, 0),
        legend.text=element_text(size=5, margin=margin(0, 4, 0, 0)),
        legend.key.size=unit(0.2, "cm"),
        legend.spacing.x=unit(0.1, "cm"),
        legend.title=element_blank())

#Plot bargraph for number of CSEPs
gg.category.CSEPs <- ggplot(genes.df, aes(y=name, x=CSEPs, fill=category)) +
  geom_bar(stat="identity", size=0.5, width=0.6) +
  labs(x="CSEPs") +
  scale_y_discrete(limits=rev(tip.order.datedtree),
                   expand=expansion(0, 0)) +
  coord_cartesian(clip="off",
                  ylim=c(0, 62)) +
  scale_x_continuous(expand=c(0, 0),
                     position="top") +
  scale_fill_manual(values=c("lightgrey", "darkgrey", "dimgrey"),
                    breaks=c("core", "accessory", "specific"),
                    labels=c("Core", "Accessory", "Specific"),
                    guide=FALSE) +
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.title.x.top=element_text(size=6, margin=margin(0, 0, 4, 0)),
        axis.text.x.top=element_text(size=3),
        panel.grid.major.x=element_line(colour="white"),
        panel.grid.minor.x=element_line(colour="white"),
        panel.grid.major.y=element_blank(),
        legend.position=c(0, 0),
        legend.text=element_text(size=5),
        legend.key.size=unit(0.2, "cm"),
        legend.spacing.x=unit(0.1, "cm"),
        legend.title=element_blank())


## REPORTED LIFESTYLES ##

#Make dataframe of reported lifestyles
lifestyle.df <- metadata[c((which(colnames(metadata) == "lifestyle") + 1):ncol(metadata))]
rownames(lifestyle.df) <- metadata$name

#Order colours
col.df$lifestyle <- factor(col.df$lifestyle,
                           levels=c("plant associate", "endophyte", "plant pathogen", "saprotroph", 
                                    "animal associate", "insect mutualist",  "animal pathogen", "human pathogen", "mycoparasite"))

#Plot dated tree
gg.maintree <- gg.datedtree.independent +
  new_scale_fill() +
  geom_nodepoint(aes(subset=(node %in% c(64, 66))),
                 size=1, colour="black") +
  geom_text_repel(data=gg.tree.data.dated.independent[gg.tree.data.dated.independent$node %in% c(64, 66),],
                  aes(x=x-max(gg.tree.data.dated.independent$x), y=y,
                      label=c("Generic concept\nsensu lato",
                              "Generic concept\nsensu stricto")),
                  hjust=1,
                  size=1.5,
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
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(5, "pt"),
        legend.text=element_text(size=4, margin=margin(0, 5, 0, 0)),
        legend.title=element_text(face="bold", size=5, margin=margin(0, 5, 0, 0)))


## GENE COPY NUMBER ##

#Filter presence-absence matrices for orthogroups present in at least one ingroup taxon
copynum.orthogroups <- orthogroups.copies.ingroup1[which(rowSums(orthogroups.copies.ingroup1) > 0),]
#copynum.CSEPs <- CSEP.count.ingroup1[which(rowSums(CSEP.count.ingroup1) > 0),]
copynum.CSEPs <- CSEP.orthogroups

#Make dataframes of gene copy number for each taxon
copynum.df.orthogroups <- data.frame(taxon=colnames(copynum.orthogroups),
                                     lifestyle=metadata$lifestyle[match(colnames(copynum.orthogroups),
                                                                        metadata$file2)],
                                     as.data.frame(t(copynum.orthogroups)))

copynum.df.CSEPs <- data.frame(taxon=colnames(copynum.CSEPs),
                                   lifestyle=metadata$lifestyle[match(colnames(copynum.CSEPs),
                                                                      metadata$file2)],
                                   as.data.frame(t(copynum.CSEPs)))

#Combine dataframes for plotting
copynum.df <- rbind(data.frame(melt(copynum.df.orthogroups), data="Orthogroups"),
                    data.frame(melt(copynum.df.CSEPs), data="CSEPs"))
#Remove copy number of 0
copynum.df <- copynum.df[copynum.df$value != 0,]

#Test for significant difference of means
tukey.copynum.orthogroups <- TukeyHSD(aov(lm(value ~ lifestyle,
                                             data=copynum.df[copynum.df$data == "Orthogroups",])))
tukey.copynum.CSEPs <- TukeyHSD(aov(lm(value ~ lifestyle,
                                           data=copynum.df[copynum.df$data == "CSEPs",])))

#Make dataframe for labelling plot
copynumlabels.df <- rbind(data.frame(tukey=multcompLetters(tukey.copynum.CSEPs[["lifestyle"]][,4])$Letters,
                                     data="CSEPs",
                                     lifestyle=names(multcompLetters(tukey.copynum.CSEPs[["lifestyle"]][,4])$Letters)),
                          data.frame(tukey=multcompLetters(tukey.copynum.orthogroups[["lifestyle"]][,4])$Letters,
                                     data="Orthogroups",
                                     lifestyle=names(multcompLetters(tukey.copynum.orthogroups[["lifestyle"]][,4])$Letters)))

#Add sample sizes for lifestyles
copynumlabels.df$num <- NA
for (i in na.omit(unique(metadata$lifestyle))) {
  copynumlabels.df$num[copynumlabels.df$lifestyle == i] <- paste0("n=", table(metadata$lifestyle[metadata$ingroup == 1])[names(table(metadata$lifestyle[metadata$ingroup == 1])) == i])
}

#Replace spaces with linebreaks
copynum.df$lifestyle <- sub(" ", "\n", copynum.df$lifestyle)
copynumlabels.df$lifestyle <- sub(" ", "\n", copynumlabels.df$lifestyle)

#Set seed for reproducible jitter
set.seed(1)

gg.copynum <- ggplot(copynum.df, aes(x=lifestyle, y=value)) +
  facet_wrap(. ~ data, scales="free",
             labeller=labeller(data=c(CSEPs="CSEPs", Orthogroups="All genes"))) +
  geom_point(position="jitter", aes(colour=lifestyle), size=0.3) +
  geom_text(data=copynumlabels.df,
            aes(x=lifestyle, y=Inf, label=tukey),
            family="mono",
            hjust=0.5,
            vjust=2,
            size=2,
            inherit.aes=FALSE) +
  geom_text(data=copynumlabels.df,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=num),
            hjust=0.5,
            vjust=2,
            size=1.5,
            show.legend=FALSE,
            inherit.aes=FALSE) +
  labs(y="Copy-number") +
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
        plot.margin=margin(2, 0, 5, 5, unit="pt"))


## STRAIN SPECIFIC GENES ##

#Filter presence-absence matrices for strain specific orthogroups
orthogroups.count.specific <- orthogroups.copies.ingroup1[which(orthogroups.stats.ingroup1$category == "specific"),]
CSEP.count.specific <- orthogroups.copies.ingroup1[
  intersect(which(orthogroups.stats.ingroup1$category == "specific"),
            which(!is.na(orthogroups.stats.ingroup1$CSEP))),]

#Make dataframe for number of strain specific orthogroups per taxon
specific.df <- data.frame(taxon=colnames(CSEP.count.specific),
                          CSEPs=NA,
                          Orthogroups=NA)
for (i in specific.df$taxon) {
  specific.df$CSEPs[specific.df$taxon == i] <- length(which(CSEP.count.specific[,i] > 0))
  specific.df$Orthogroups[specific.df$taxon == i] <- length(which(orthogroups.count.specific[,i] > 0))
}
#Add lifestyle
specific.df$lifestyle <- metadata$lifestyle[match(specific.df$taxon, metadata$file2)]

#Test for significant difference in means
tukey.CSEPs <- TukeyHSD(aov(lm(CSEPs ~ lifestyle, data=specific.df)))
tukey.orthogroups <- TukeyHSD(aov(lm(Orthogroups ~ lifestyle, data=specific.df)))
#Make dataframe for labelling plot
specificlabels.df <- rbind(data.frame(tukey=multcompLetters(tukey.CSEPs[["lifestyle"]][,4])$Letters,
                                      variable="CSEPs",
                                      lifestyle=names(multcompLetters(tukey.CSEPs[["lifestyle"]][,4])$Letters)),
                           data.frame(tukey=multcompLetters(tukey.orthogroups[["lifestyle"]][,4])$Letters,
                                      variable="Orthogroups",
                                      lifestyle=names(multcompLetters(tukey.orthogroups[["lifestyle"]][,4])$Letters)))
#Melt dataframe for plotting
specific.df <- melt(specific.df)
#Replace spaces with linebreaks
specific.df$lifestyle <- sub(" ", "\n", specific.df$lifestyle)
specificlabels.df$lifestyle <- sub(" ", "\n", specificlabels.df$lifestyle)
#Add sample size
specificlabels.df$num <- copynumlabels.df$num[match(specificlabels.df$lifestyle, copynumlabels.df$lifestyle)]

#Plot boxplot of strain specific orthogroups across lifestyles
gg.specific <- ggplot(specific.df, aes(x=lifestyle, y=value, fill=lifestyle)) +
  facet_wrap(. ~ variable, scales="free",
             labeller=labeller(variable=c(CSEPs="CSEPs", Orthogroups="All genes"))) +
  geom_violin(fill="white", colour="grey", lty="dotted", size=0.3) +
  geom_boxplot(width=0.2, size=0.3, outlier.size=0.5) +
  geom_text(data=specificlabels.df,
            aes(x=lifestyle, y=Inf, label=tukey),
            family="mono",
            hjust=0.5,
            vjust=2,
            size=2,
            inherit.aes=FALSE) +
  geom_text(data=specificlabels.df,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=num),
            hjust=0.5,
            vjust=2,
            size=1.5,
            inherit.aes=FALSE) +
  labs(y="Strain specific genes") +
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
  guides(fill=guide_legend(direction="horizontal",
                           nrow=2),
         colour="none") +
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
        plot.margin=margin(0, 0, 5, 0, unit="pt"))

#Write boxplot of strain specific genes to file (Supplementary Figure 5)
tiff(file=paste0("SupplementaryFig5-", Sys.Date(), ".tiff"),
     height=2, width=3.4, unit="in", res=600, compression="lzw")
gg.specific
dev.off()


## LIFESTYLE COMPARISON ##

#For CSEPs and all orthogroups...
for (i in c("CSEPs", "orthogroups")){
  
  print(i)
  #Read in lifestyle test results
  phy.pca.result <- read.csv(paste0("lifestyle_comparison/", i, "/metadata.csv"))
  lifestyle.data <- read.csv(paste0("lifestyle_comparison/", i, "/data.csv"), row.names="genome")
  #Make distance matrix of orthogroup content
  dist <- vegdist(lifestyle.data, method="jaccard")
  #Do permanova
  permanova <- adonis2(formula=dist ~ PC1 + PC2 + lifestyle, data=phy.pca.result, permutations=9999)
  
  print(paste0("Phylogeny: ", round(sum(permanova$R2[1:2]) * 100), "%"))
  print(paste0("Lifestyle: ", round(sum(permanova$R2[3]) * 100), "%"))
  
  assign(paste0("permanova.", i), permanova)
  
}

#Make dataframe with pairwise PERMANOVA results
pw.lifestyle.genes <- rbind(data.frame(melt(as.matrix(read.csv(
  "lifestyle_comparison/CSEPs/pairwiseComparisons.csv",
  row.names=1)), na.rm=TRUE), data="CSEPs"),
  data.frame(melt(as.matrix(read.csv(
    "lifestyle_comparison/orthogroups/pairwiseComparisons.csv",
    row.names=1)), na.rm=TRUE), data="Orthogroups"))

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
                                        "%\nLifestyle: ", round(sum(permanova.CSEPs$R2[3]) * 100), "%"),
                                 paste0("Phylogeny: ", round(sum(permanova.orthogroups$R2[1:2]) * 100),
                                        "%\nLifestyle: ", round(sum(permanova.orthogroups$R2[3]) * 100), "%")),
                           data=c("CSEPs", "Orthogroups"))


#Plot grid
gg.pwperm <- ggplot(pw.lifestyle.genes, aes(Var2, Var1, fill=value>0.05)) +
  facet_grid(. ~ data, labeller=labeller(data=c(CSEPs="CSEPs", Orthogroups="All genes"))) +
  geom_tile(color="grey", size=1, alpha=0.5, show.legend=FALSE) +
  geom_text(aes(label=label, colour=value>0.05), size=1.5, show.legend=FALSE) +
  annotate("text", x=4, y=2, label="PERMANOVA", size=1.5, fontface="bold") +
  geom_text(data=permanova.df, aes(x=4, y=1.5, label=lab), size=1.5, inherit.aes=FALSE) +
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
        plot.margin=unit(c(0, 5, 0, 0), "pt")) +
  coord_fixed()

## COMBINE ##

#Write to file
tiff(file=paste0("Fig3-", Sys.Date(), ".tiff"),
     height=6, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.lifestyles.grid, gg.category.CSEPs, gg.category.orthogroups,
                    nrow=1, rel_widths=c(6, 1, 1), align="h", axis="bt"),
          plot_grid(gg.pwperm, gg.copynum, nrow=1, rel_widths=c(1, 1), labels=c("", "C"), label_size=10),
          nrow=2,
          rel_heights=c(2, 1),
          labels="AUTO",
          label_size=10)
dev.off()


############################    FIGURE 4    ####################################
############################   SELECTION    ####################################

## BUSTED AND aBSREL ##

#Vector of core, single-copy orthogroups
core.SC.orthogroups <- Reduce(intersect,
                              list(orthogroups.stats.ingroup0$orthogroup[which(
                                orthogroups.stats.ingroup0$copy_number == "single")],
                                orthogroups.stats.ingroup0$orthogroup[which(
                                  orthogroups.stats.ingroup0$category == "core")]))

#Core, single-copy CSEPs
core.SC.mixed <- Reduce(intersect,
                        list(orthogroups.stats.ingroup0$orthogroup[which(
                          orthogroups.stats.ingroup0$copy_number == "single")],
                          orthogroups.stats.ingroup0$orthogroup[which(
                            orthogroups.stats.ingroup0$category == "core")],
                          orthogroups.stats.ingroup0$orthogroup[which(
                            orthogroups.stats.ingroup0$CSEP != "")]))

#Make dataframe for selection results
selection.df <- data.frame(orthogroup=core.SC.orthogroups,
                           CSEP="N", busted="N", absrel="N", busted.success="Y", absrel.success="Y")
#Orthogroups predicted to be CSEPs
selection.df$CSEP[match(core.SC.mixed, selection.df$orthogroup)] <- "Y"

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
      absrel.p[[names(absrel.results[["branch attributes"]][["0"]][j])]][i] <- absrel.results[["branch attributes"]][["0"]][[names(absrel.results[["branch attributes"]][["0"]])[j]]][["Corrected P-value"]]
    }
    
  } else {
    selection.df$absrel.success[match(i, selection.df$orthogroup)] <- "N"
  }
  
}

#Filter for significant aBSREL results
absrel.p.sig <- lapply(absrel.p, function(x) names(which(x < 0.05)))

#Add aBSREL results to selection dataframe 
selection.df$absrel[which(selection.df$orthogroup %in% unlist(absrel.p.sig))] <- "Y"

#Make matrix for euler plot of selection results
selection.mat.labels <- c("Total"=length(union(which(selection.df$busted.success == "Y"),
                                               which(selection.df$absrel.success == "Y"))),
                          "BUSTED"=0,
                          "aBSREL"=0,
                          "CSEPs"=0,
                          "BUSTED&aBSREL"=0,
                          "CSEPs&BUSTED"=0,
                          "CSEPs&aBSREL"=0,
                          "Total&BUSTED"=length(Reduce(intersect,
                                                       list(which(selection.df$CSEP == "N"),
                                                            which(selection.df$busted == "Y"),
                                                            which(selection.df$absrel == "N")))),
                          "Total&aBSREL"=length(Reduce(intersect,
                                                       list(which(selection.df$CSEP == "N"),
                                                            which(selection.df$busted == "N"),
                                                            which(selection.df$absrel == "Y")))),
                          "Total&CSEPs"=length(Reduce(intersect,
                                                      list(which(selection.df$CSEP == "Y"),
                                                           which(selection.df$busted == "N"),
                                                           which(selection.df$absrel == "N")))),
                          "Total&BUSTED&aBSREL"=length(Reduce(intersect,
                                                              list(which(selection.df$CSEP == "N"),
                                                                   which(selection.df$busted == "Y"),
                                                                   which(selection.df$absrel == "Y")))),
                          "Total&BUSTED&CSEPs"=length(Reduce(intersect,
                                                             list(which(selection.df$CSEP == "Y"),
                                                                  which(selection.df$busted == "Y"),
                                                                  which(selection.df$absrel == "N")))),
                          "Total&aBSREL&CSEPs"=length(Reduce(intersect,
                                                             list(which(selection.df$CSEP == "Y"),
                                                                  which(selection.df$busted == "N"),
                                                                  which(selection.df$absrel == "Y")))),
                          "Total&BUSTED&aBSREL&CSEPs"=length(Reduce(intersect,
                                                                    list(which(selection.df$CSEP == "Y"),
                                                                         which(selection.df$busted == "Y"),
                                                                         which(selection.df$absrel == "Y")))))

#Number of core SC orthogroups positively selected (consensus between aBSREL and BUSTED)
print(paste0(selection.mat.labels[["Total&BUSTED&aBSREL&CSEPs"]] + selection.mat.labels[["Total&BUSTED&aBSREL"]],
      ", ",
      round((selection.mat.labels[["Total&BUSTED&aBSREL&CSEPs"]] + selection.mat.labels[["Total&BUSTED&aBSREL"]]) /
        selection.mat.labels[["Total"]] * 100),
      "%"))

#Make second matrix with altered numbers to fix plot
selection.mat <- selection.mat.labels
selection.mat[["Total&BUSTED&CSEPs"]] <-selection.mat[["Total&aBSREL&CSEPs"]]
selection.mat[["Total&CSEPs"]] <- selection.mat[["Total&BUSTED&aBSREL&CSEPs"]]

#Create and plot euler diagram
set.seed(1)
selection.euler <- euler(selection.mat, shape="ellipse")
eulerr_options(padding=unit(0.5, "pt"))
euler <- plot(selection.euler,
              fills=list(fill=c("white", "dimgrey", "dimgrey", "white", rep("", 6),
                                "#56B4E9", "grey", "grey", rep("", 1), "#56B4E9"),
                         alpha=0.3),
              labels=list(labels=c("Genes",
                                   "BUSTED",
                                   "aBSREL",
                                   "CSEPs"),
                          cex=0.35,
                          col="dimgrey",
                          padding=unit(c(1, 1), "pt")),
              edges=list(col="dimgrey", lty=c("solid", "dotted", "dotted", "solid"), lwd=0.5),
              quantities=list(labels=c(selection.mat.labels[["Total"]],
                                       rep("", 3),
                                       selection.mat.labels[["Total&BUSTED"]],
                                       selection.mat.labels[["Total&aBSREL"]],
                                       selection.mat.labels[["Total&CSEPs"]], rep("", 3),
                                       paste0("Positively\nselected\nconsensus\n\n\n",
                                              selection.mat.labels[["Total&BUSTED&aBSREL"]]), 
                                       selection.mat.labels[["Total&BUSTED&CSEPs"]],
                                       selection.mat.labels[["Total&aBSREL&CSEPs"]], "",
                                       selection.mat.labels[["Total&BUSTED&aBSREL&CSEPs"]]),
                              col=c(rep("dimgrey", 4), "black", rep("dimgrey", 2), "black"),
                              cex=0.3))

#Convert to grob
euler.grob <- as.grob(euler)

#Read in dataframe matching aBSREL tree nodes to IQ-TREE species tree nodes
absrel.df <- read.csv("selection/hyphy/absrel_nodes.csv")
#Add tip labels
absrel.df <- rbind(absrel.df, data.frame(node=1:length(iqtree$tip.label), absrel=iqtree$tip.label))

##Add number of significant aBSREL tests per taxon that are also supported by busted
absrel.p.count <- lengths(lapply(absrel.p.sig, function(x)
  x[which(x %in% selection.df$orthogroup[which(selection.df$busted == "Y")])]))
absrel.df$num <- absrel.p.count[match(absrel.df$absrel, names(absrel.p.count))]

#Add names of positively selected CSEPs
absrel.p.names <- unlist(lapply(absrel.p.sig, function(x)
  paste(x[intersect(which(x %in% core.SC.mixed),
                    which(x %in% selection.df$orthogroup[which(selection.df$busted == "Y")]))],
          collapse=" ")))
absrel.p.names <- absrel.p.names[absrel.p.names != ""]
absrel.p.names <- gsub("OG000", "", absrel.p.names)

absrel.df$CSEPs <- NA

absrel.df$CSEPs <- absrel.p.names[match(absrel.df$absrel, names(absrel.p.names))]
absrel.df$CSEPs[absrel.df$CSEPs == ""] <- NA
absrel.df$CSEPs <- gsub('(?=(?:.{10})+$)', "\n", absrel.df$CSEPs, perl = TRUE)

#Plot species tree with no branch lengths
gg.selection <- ggtree(iqtree.tree, branch.length="none") %<+% metadata

#Create vector with the order of the tree tips for plotting
tip.order.speciestree <- rev(gg.tree.data.iqtree$label[gg.tree.data.iqtree$isTip == "TRUE"])
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
tiplabel.face.speciestree <- rep("italic", length(iqtree.tree$tip.label))
tiplabel.face.speciestree[which(metadata$new[match(iqtree.tree$tip.label, metadata$name)] == "Y")] <- "bold.italic"

#Make dataframe of species complex nodes
sc.df.iq <- data.frame(sc=names(table(metadata$speciescomplex.abb[match(iqtree.tree$tip.label, metadata$name)]))[table(metadata$speciescomplex.abb[match(iqtree.tree$tip.label, metadata$name)]) > 1],
                       node=NA)

#Get nodes for each species complex
for (i in 1:length(sc.df.iq$sc)) {
  sc.df.iq$node[i] <- getMRCA(iqtree.tree, metadata$name[metadata$speciescomplex.abb == sc.df.iq$sc[i]])
}

sc.df.iq2 <- data.frame(sc=unique(metadata$speciescomplex.abb)[is.na(match(unique(metadata$speciescomplex.abb),
                                                                           sc.df.iq$sc))],
                     node=NA)

sc.df.iq2$node <- gg.tree.data.iqtree$node[match(sc.df.iq2$sc, gg.tree.data.iqtree$speciescomplex.abb)]

sc.df.iq <- rbind(sc.df.iq, sc.df.iq2)
sc.df.iq <- sc.df.iq[match(na.omit(unique(gg.tree.data.iqtree$speciescomplex.abb)), sc.df.iq$sc),]
sc.df.iq$box <- rep(c(0,1), length.out=length(sc.df.iq$sc))

for (i in 1:length(sc.df.iq$sc)) {
  
  if (sc.df.iq$box[i] == 1) {
    
    gg.selection <- gg.selection +
      geom_highlight(node=sc.df.iq$node[sc.df.iq$sc == sc.df.iq$sc[i]], fill="#000000", alpha=0.075, extend=10)
    
  } else {
    
    gg.selection <- gg.selection +
      geom_highlight(node=sc.df.iq$node[sc.df.iq$sc == sc.df.iq$sc[i]], fill="#000000", alpha=0.04, extend=10)
    
  }
  
}

#Write numbered selection tree to file (Supplementary Figure 6)
tiff(file=paste0("SupplementaryFig6-", Sys.Date(), ".tiff"),
     height=4, width=6.75, unit="in", res=600, compression="lzw")

gg.selection %<+% absrel.df +
  geom_tree(size=1, aes(colour=num)) +
  xlim(0, 22) +
  scale_colour_gradient2(trans="pseudo_log",
                         low="black", mid="#F0E442", high="#CC79A7", midpoint=1,
                         breaks=c(0, 10, 100),
                         limits=c(0, 100),
                         labels=c(0 , 10, 100),
                         na.value="dimgrey",
                         guide=guide_colourbar(title="Genes showing\npositive selection",
                                               title.position="top",
                                               direction="horizontal")) +
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
  geom_tiplab(aes(label=tiplab), fontface=tiplabel.face.speciestree, size=1.5, offset=0.05, colour="black") +
  coord_cartesian(clip="off") +
  annotation_custom(euler.grob, xmin=0, xmax=7, ymin=40, ymax=62) +
  theme(legend.box="vertical",
        legend.title=element_text(size=5, face="bold"),
        legend.text=element_text(size=5, margin=margin(0, 3, 0, 0)),
        legend.position=c(0.13, 0.45),
        legend.key.size=unit(6, "pt"),
        legend.spacing=unit(0, "pt"),
        plot.margin=margin(0, 0, 5, -16, unit="pt"))

dev.off()

#Plot aBSREL results on tree
gg.selection2 <- gg.selection %<+% absrel.df +
  geom_tree(size=1, aes(colour=num)) +
  xlim(0, 22) +
  scale_colour_gradient2(trans="pseudo_log",
                         low="black", mid="#F0E442", high="#CC79A7", midpoint=1,
                         breaks=c(0, 10, 100),
                         limits=c(0, 100),
                         labels=c(0 , 10, 100),
                         na.value="dimgrey",
                         guide=guide_colourbar(title="Genes showing\npositive selection",
                                               title.position="top",
                                               direction="horizontal")) +
  geom_label2(aes(x=branch, label=CSEPs, fill=num),
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
  geom_tiplab(aes(label=tiplab), fontface=tiplabel.face.speciestree, size=1.5, offset=0.05, colour="black") +
  new_scale_colour() +
  geom_tippoint(aes(colour=lifestyle), size=0.8) +
  scale_colour_manual(values=col.df$colour[na.omit(match(sort(unique(metadata$lifestyle[match(iqtree.tree$tip.label, metadata$name)])), col.df$lifestyle))],
                      labels=str_to_sentence(col.df$lifestyle[na.omit(match(sort(unique(metadata$lifestyle[match(iqtree.tree$tip.label, metadata$name)])), col.df$lifestyle))]),
                      na.translate=FALSE,
                      guide=guide_legend(title="Lifestyle",
                                         ncol=2,
                                         title.hjust=0,
                                         order=1)) +
  coord_cartesian(clip="off") +
  annotation_custom(euler.grob, xmin=0, xmax=7, ymin=40, ymax=62) +
  theme(legend.box="vertical",
        legend.title=element_text(size=5, face="bold"),
        legend.text=element_text(size=5, margin=margin(0, 3, 0, 0)),
        legend.position=c(0.13, 0.5),
        legend.key.size=unit(6, "pt"),
        legend.spacing=unit(0, "pt"),
        plot.margin=margin(0, 0, 5, -16, unit="pt"))


## CONTRAST-FEL ##

#Make dataframe for Contrast-FEL results
contrastfel.df <- data.frame(ortho=rep(core.SC.orthogroups, each=length(na.omit(unique(metadata$lifestyle)))),
                             lifestyle=rep(na.omit(unique(metadata$lifestyle)), times=length(core.SC.orthogroups)),
                             increase=NA,
                             decrease=NA)

contrastfel.df$CSEP <- orthogroups.stats.ingroup1$CSEP[match(contrastfel.df$ortho,
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

#Test for significant difference in means
tukey.increase <- TukeyHSD(aov(lm(increase ~ lifestyle, data=contrastfel.df[which(contrastfel.df$increase > 0),])))
tukey.decrease <- TukeyHSD(aov(lm(decrease ~ lifestyle, data=contrastfel.df[which(contrastfel.df$decrease > 0),])))
#Make dataframe for labelling plot
sitelabels.df <- rbind(data.frame(tukey=multcompLetters(tukey.increase[["lifestyle"]][,4])$Letters,
                                  variable="increase",
                                  lifestyle=names(multcompLetters(tukey.increase[["lifestyle"]][,4])$Letters)),
                       data.frame(tukey=multcompLetters(tukey.decrease[["lifestyle"]][,4])$Letters,
                                  variable="decrease",
                                  lifestyle=names(multcompLetters(tukey.decrease[["lifestyle"]][,4])$Letters)))

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

#Melt dataframe for plotting
contrastfel.sites <- melt(contrastfel.df)
#Remove genes with no sites
contrastfel.sites <- contrastfel.sites[which(contrastfel.sites$value > 0),]

#Add column for whether the positive selection according to aBSREL/BUSTED
contrastfel.sites$positive.selection <- NA
contrastfel.sites$positive.selection[
  intersect(which(!is.na(contrastfel.sites$CSEP)),
            which(contrastfel.sites$value != 0))][
              contrastfel.sites$ortho[intersect(which(!is.na(contrastfel.sites$CSEP)),
                                                which(contrastfel.sites$value != 0))] %in%
                selection.df$orthogroup[intersect(which(selection.df$busted == "Y"),
                                                  which(selection.df$absrel == "Y"))]] <- "Y"

#Remove prefix
contrastfel.sites$ortho <- sub("OG000", "", contrastfel.sites$ortho)
#Replace spaces with linebreaks
contrastfel.sites$lifestyle <- sub(" ", "\n", contrastfel.sites$lifestyle)

#Plot violin plot of Contrast-FEL results
gg.siterates <- ggplot(contrastfel.sites, aes(x=lifestyle, y=value, fill=lifestyle)) +
  facet_wrap(. ~ variable, labeller=labeller(variable=c(increase="Higher relative selective pressure",
                                                        decrease="Lower relative selective pressure"))) +
  geom_violin(aes(colour=lifestyle), size=0.3) +
  geom_point(data=contrastfel.sites[which(contrastfel.sites$positive.selection == "Y"),],
             size=0.2,
             colour="black") +
  geom_text(data=sitelabels.df,
            aes(x=lifestyle, y=Inf, label=tukey),
            family="mono",
            hjust=0.5,
            vjust=2,
            size=2,
            inherit.aes=FALSE) +
  geom_text(data=sitelabels.df,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=paste0("n=", num)),
            hjust=0.5,
            vjust=6,
            size=1.5,
            inherit.aes=FALSE) +
  geom_text_repel(data=contrastfel.sites[which(contrastfel.sites$positive.selection == "Y"),],
                  aes(label=ortho),
                  segment.size=0.2,
                  min.segment.length=unit(0, 'lines'),
                  segment.colour="grey",
                  colour="black",
                  nudge_x=0.5,
                  direction="y",
                  max.overlaps=Inf,
                  size=1) +
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
        plot.margin=margin(10, 0, 10, 10, unit="pt"))

## COMBINE ##

#Write to file
tiff(file=paste0("Fig4-", Sys.Date(), ".tiff"),
     height=6, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(gg.selection2,
          gg.siterates,
          nrow=2,
          rel_heights=c(2, 1),
          labels="AUTO",
          label_size=10)
dev.off()


############################    FIGURE 5    ####################################
########################   CODON OPTIMISATION    ###############################

#Load codon optimisation results
load("selection/codon_optimisation/codon_optimisation-2021-09-10.RData")

## LIFESTYLE BOXPLOT ##

#Add lifestyle and name to optimisation (S values) dataframe
s.df$lifestyle <- metadata$lifestyle[match(s.df$taxon, metadata$file2)]
s.df$name <- metadata$name[match(s.df$taxon, metadata$file2)]
s.df <- s.df %>% select(name, everything())

#Remove outgroup
s.df <- s.df[-which(is.na(s.df$lifestyle)),]

#Test for significant difference in optimisation
tukey.codons <- TukeyHSD(aov(lm(S ~ lifestyle, data=s.df)))

#Make dataframe for labelling plot
codonlabels.df <- data.frame(tukey=multcompLetters(tukey.codons[["lifestyle"]][,4])$Letters,
                             lifestyle=names(multcompLetters(tukey.codons[["lifestyle"]][,4])$Letters))

#Replace spaces with linebreaks
s.df$lifestyle <- sub(" ", "\n", s.df$lifestyle)
codonlabels.df$lifestyle <- sub(" ", "\n", codonlabels.df$lifestyle)

#Add columns for statistical test results
codonlabels.df$sig <- NA
codonlabels.df$pos <- NA

#For each lifestyle
for (i in unique(s.df$lifestyle)) {
  
  #Do Wilcox test for differences between CSEPs and non-CSEPs
  s.p <- wilcox.test(x=s.df$S.CSEP[s.df$lifestyle == i], y=s.df$S.other[s.df$lifestyle == i])$p.value
  
  print(paste(i, "=", signif(s.p, digits=1)))
  
  #If the p value is significant...
  if (s.p <= 0.05) {
    
    #Add an asterisk label
    codonlabels.df$sig[codonlabels.df$lifestyle == i] <- "*"
    #Add position for label
    codonlabels.df$pos[codonlabels.df$lifestyle == i] <-max(c(s.df$S[which(s.df$lifestyle == i)],
                                                              s.df$S.CSEP[which(s.df$lifestyle == i)],
                                                              s.df$S.other[which(s.df$lifestyle == i)]))
    
  }
  
}

#Add sample size
codonlabels.df$num <- copynumlabels.df$num[match(codonlabels.df$lifestyle, copynumlabels.df$lifestyle)]

#Melt dataframe for plotting
s.df2 <- melt(s.df[-3])

#Plot boxplot of codon optimisation
gg.codons <- ggplot(s.df2, aes(x=lifestyle, y=value, fill=lifestyle, alpha=variable)) +
  geom_boxplot(position=position_dodge(width=0.4), width=0.3, size=0.3, outlier.size=0.5) +
  geom_text(data=codonlabels.df,
            aes(x=lifestyle, y=pos, label=sig),
            fontface="bold",
            vjust=-0.1,
            size=3,
            inherit.aes=FALSE) +
  geom_text(data=codonlabels.df,
            aes(x=lifestyle, y=Inf, label=tukey),
            family="mono",
            hjust=0.5,
            vjust=2,
            size=2,
            inherit.aes=FALSE) +
  geom_text(data=codonlabels.df,
            aes(x=lifestyle, y=-Inf, colour=lifestyle, label=num),
            hjust=0.5,
            vjust=6,
            size=1.5,
            inherit.aes=FALSE) +
  labs(y="Codon optimisation (S)") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
  scale_colour_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                      labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                    col.df$lifestyle)]),
                      name=NULL,
                      guide=FALSE) +
  scale_fill_manual(values=col.df$colour[match(sort(unique(na.omit(metadata$lifestyle))), col.df$lifestyle)],
                    labels=str_to_sentence(col.df$lifestyle[match(sort(unique(na.omit(metadata$lifestyle))),
                                                                  col.df$lifestyle)]),
                    name=NULL) +
  scale_alpha_manual(values=c(1, 0.3),
                     labels=c("CSEPs", "Other"),
                     name=NULL) +
  coord_cartesian(clip="off") +
  guides(alpha=guide_legend(direction="horizontal", override.aes=list(fill="dimgrey")),
         fill=FALSE) +
  theme_minimal() +
  theme(legend.position="top",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.text=element_text(size=6, margin=margin(0, 5, 0, 0)),
        axis.text.x=element_text(colour=c("#009E73","#56B4E9", "#D55E00", "#9AE324", "dimgrey", "#0072B2"),
                                 size=4, face="bold", vjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        plot.margin=margin(0, 10, 10, 10, unit="pt"))


## RANGE DOTPLOT ##

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
range.pgls.pagel <- gls(S ~ range,
                        correlation=corPagel(1, phy=plgs.range.tree, form=~name),
                        data=range.s.df, method="ML")
range.pgls.brownian <- gls(S ~ range,
                           correlation=corBrownian(phy=plgs.range.tree, form=~name),
                           data=range.s.df, method="ML")
range.pgls.blomberg <- gls(S ~ range, 
                           correlation=corBlomberg(1, phy=plgs.range.tree, form=~name, fixed=TRUE),
                           data=range.s.df, method="ML")

#Select the best model according to AIC 
selected.model <- rownames(anova(range.pgls.pagel, range.pgls.brownian, range.pgls.blomberg))[which.min(anova(range.pgls.pagel, range.pgls.brownian, range.pgls.blomberg)$AIC)]
#Add results for selected model to dataframe
range.s.df$pgls <- predict(get(selected.model))

#Fix bug in dotplots: https://stackoverflow.com/questions/51406933/how-to-use-ggplot2s-geom-dotplot-with-symmetrically-spaced-and-separated-dots
dotstackGrob <- function(
  x = unit(0.5, "npc"),     # x pos of the dotstack's origin
  y = unit(0.5, "npc"),     # y pos of the dotstack's origin
  stackaxis = "y",
  dotdia = unit(1, "npc"),  # Dot diameter in the non-stack axis, should be in npc
  stackposition = 0,        # Position of each dot in the stack, relative to origin
  stackratio = 1,           # Stacking height of dots (.75 means 25% dot overlap)
  default.units = "npc", name = NULL, gp = gpar(), vp = NULL)
{
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  if (!is.unit(dotdia))
    dotdia <- unit(dotdia, default.units)
  if (attr(dotdia,"unit") != "npc")
    warning("Unit type of dotdia should be 'npc'")
  
  grob(x = x, y = y, stackaxis = stackaxis, dotdia = dotdia,
       stackposition = stackposition, stackratio = stackratio,
       name = name, gp = gp, vp = vp, cl = "dotstackGrob")
}

#' @export
makeContext.dotstackGrob <- function(x, recording = TRUE) {
  # Need absolute coordinates because when using npc coords with circleGrob,
  # the radius is in the _smaller_ of the two axes. We need the radius
  # to instead be defined in terms of the non-stack axis.
  xmm <- convertX(x$x, "mm", valueOnly = TRUE)
  ymm <- convertY(x$y, "mm", valueOnly = TRUE)
  
  if (x$stackaxis == "x") {
    dotdiamm <- convertY(x$dotdia, "mm", valueOnly = TRUE)
    xpos <- xmm + dotdiamm * (x$stackposition * x$stackratio)
    ypos <- ymm
  } else if (x$stackaxis == "y") {
    dotdiamm <- convertX(x$dotdia, "mm", valueOnly = TRUE)
    xpos <- xmm
    ypos <- ymm + dotdiamm * (x$stackposition * x$stackratio + (1 - x$stackratio) / 2)
  }
  
  circleGrob(
    x = xpos, y = ypos, r = dotdiamm / 2, default.units = "mm",
    name = x$name, gp = x$gp, vp = x$vp
  )
}

#Plot dotplot of range against codon optimisation
gg.range <- ggplot(range.s.df, aes(x=range, y=S)) +
  geom_dotplot(aes(x=as.factor(range)),
               stackratio=2.5,
               binaxis="y", stackdir="center", dotsize=0.5, fill="dimgrey", colour="dimgrey") +
  geom_line(aes(y=pgls)) +
  geom_smooth(method="lm", colour="black", linetype="dashed", size=0.5, se=FALSE) +
  geom_text(data=rangelabels.df,
            aes(x=range, y=-Inf, label=paste0("n=", n)),
            colour="dimgrey",
            hjust=0.5,
            vjust=5,
            size=1.5,
            inherit.aes=FALSE) +
  labs(x="Number of reported lifestyles", y="Codon optimisation (S)") +
  annotate("text", x=4, y=0.57, label="PGLS", fontface="bold", size=2) +
  annotate("text", x=4, y=0.55, label=paste0("p=", signif(anova(range.pgls.brownian)[2, 3], 1)), size=2) +
  annotate("text", x=5.5, y=0.57, label="Pearson's\ncorrelation", fontface="bold", size=2) +
  annotate("text", x=5.5, y=0.53, label=paste0("r=", range.corr.r, "\np=", range.corr.p), size=2) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(size=5, face="bold"),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=6, margin=margin(12, 0, 0, 0)),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=4),
        plot.margin=margin(10, 10, 0, 10, unit="pt"))

#Write lifestyle and range boxplots to file (Supplementary Figure 7)
tiff(file=paste0("SupplementaryFig7-", Sys.Date(), ".tiff"),
     height=4, width=3.4, unit="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.codons, gg.range, ncol=1, labels="AUTO", label_size=10))
dev.off()


## CORRELATION WITH PHYLOGENY ##

#Read in and format phylogenetic distance matrix from lifestyle comparison test
phyl.dist <- read.csv("lifestyle_comparison/orthogroups/phyldistmatrix.csv")
rownames(phyl.dist) <- phyl.dist$X
phyl.dist$X <- NULL

#Recreate PCA of phylogenetic distances
phylpca <- prcomp(phyl.dist, rank=2)
pca.df <- as.data.frame(phylpca$x)
#Add species complex and S values
pca.df$speciescomplex.abb <- metadata$speciescomplex.abb[match(rownames(pca.df), metadata$short.tip)]
pca.df$s.other <- s.df$S.other[match(metadata$name[match(rownames(pca.df), metadata$short.tip)], s.df$name)]
pca.df$s.CSEP <- s.df$S.CSEP[match(metadata$name[match(rownames(pca.df), metadata$short.tip)], s.df$name)]

#For CSEPs and non-CSEPs...
for (i in c("CSEP", "other")) {
  
  #Fit S values to PCA
  pca.ordi.s <- ordisurf(as.formula(paste0("phylpca ~ s.", i)), data=pca.df, plot=FALSE)
  #Pull out coordinates for plotting
  pca.ordi.s.gg <- expand.grid(x=pca.ordi.s$grid$x, y=pca.ordi.s$grid$y)
  #Get z scores
  pca.ordi.s.gg$z <- as.vector(pca.ordi.s$grid$z)
  #Remova NAs
  pca.ordi.s.gg <- data.frame(na.omit(pca.ordi.s.gg))
  
  assign(paste0("pca.ordi.s.", i), pca.ordi.s)
  assign(paste0("pca.ordi.s.gg.", i), pca.ordi.s.gg)
  
}

#Combine S coordinates into one dataframe for plotting
pca.ordi.df <- rbind(data.frame(data="CSEPs",
                                pca.ordi.s.gg.CSEP),
                     data.frame(data="Other",
                                pca.ordi.s.gg.other))

#Calculate centroids for each species complex
pca.centroids.df <- aggregate(. ~ speciescomplex.abb, pca.df[1:3], mean)

#Make dataframe for labelling plot
pcalabels.df <- data.frame(R=c(signif(summary(pca.ordi.s.CSEP)$r.sq, 1),
                               signif(summary(pca.ordi.s.other)$r.sq, 1)),
                           p=c(signif(summary(pca.ordi.s.CSEP)$s.table[4], 1),
                               signif(summary(pca.ordi.s.other)$s.table[4], 1)),
                           data=c("CSEPs", "Other"))

#Plot PCA of phylogenetic distances with S values fitted
gg.pca <- ggplot(pca.centroids.df, aes(x=PC1, y=PC2,
                             colour=speciescomplex.abb, fill=speciescomplex.abb, shape=speciescomplex.abb)) +
  facet_wrap(~ data, scales="free", ncol=1) +
  geom_contour(data=pca.ordi.df, 
               aes(x=x, y=y, z=z),
               colour="dimgrey",
               show.legend=FALSE,
               size=0.3,
               inherit.aes = FALSE) +
  geom_label_contour(data=pca.ordi.df,
                     aes(x=x, y=y, z=z),
                     colour="dimgrey",
                     binwidth=0.05,
                     size=2,
                     label.size=NA,
                     label.padding=unit(0.1, "lines"),
                     inherit.aes = FALSE) +
  geom_point(size=1.5, stroke=0.5, position=position_jitter(width=0.05, height=0.05, seed=1)) +
  scale_shape_manual(values=c(1:13)) +
  annotate("text", x=-1, y=-1.3, label="ordisurf", fontface="bold", size=1.5) +
  geom_text(data=pcalabels.df,
            aes(x=-1, y=-1.6, 
                label=paste0("adj. R=", R,
                             "\np=", p)),
            size=1.5,
            inherit.aes=FALSE) +
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
s.df$genus <- "Fusarium"
s.df$genus[match(metadata$name[grep("'", metadata$speciescomplex.abb)], s.df$name)] <- "Allied"

#Test for significant difference in overall S value between genera
t.test(s.df$S[s.df$genus == "Allied"], s.df$S[s.df$genus == "Fusarium"])

#Plot boxplot of S values between genera
gg.codongenus <- ggplot(s.df, aes(x=genus, y=S)) +
  geom_violin(colour="grey", lty="dashed", size=0.3) +
  geom_boxplot(width=0.2, size=0.3, outlier.size=0.3) +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text.y=element_text(size=2.5, margin=margin(0, 0, 0, 0)),
        axis.text.x=element_text(size=3.5, face=c("plain", "italic"), margin=margin(0, 0, 0, 0)),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        plot.margin=margin(3, 3, 3, 3),
        plot.background=element_rect(colour="black", fill="white", size=0.7))

#Remove outgroup from highlighting dataframe
sc.df.pca <- sc.df.iq[!sc.df.iq$sc == "outgroup",]

#Plot tree for PCA legend
gg.pcatree <- ggtree(iqtree.tree, branch.length="none", linetype=NA) %<+% metadata +
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



#Write to file
#tiff(file=paste0("Fig5-", Sys.Date(), ".tiff"),
#     height=4, width=3.25, unit="in", res=600, compression="lzw")
#plot_grid(gg.pca, gg.pcatree, rel_widths=c(1.5, 1), align="h", axis="tb")
#dev.off()

#Write to file
#tiff(file=paste0("Fig5-", Sys.Date(), ".tiff"),
#     height=4, width=6.75, unit="in", res=600, compression="lzw")
#plot_grid(plot_grid(gg.codons, gg.range, ncol=1, labels="AUTO", label_size=10),
#          plot_grid(gg.pca,
 #                   gg.pcatree +
  #                    annotation_custom(ggplotGrob(gg.codongenus), ymin=55, ymax=65, xmin=-25, xmax=-45),
   #                 rel_widths=c(1.5, 1), align="h", axis="tb",
    #                labels=c("C", ""), label_size=10),
     #     ncol=2)
#dev.off()


## CODON USAGE BIAS HIERARCHICAL CLUSTERING ##

#Combine RSCU results for all taxa
rscu.grid <- t(cbind(as.data.frame(mget(ls(pattern="rscu.")))))
rownames(rscu.grid) <- metadata$name[match(sub("rscu.", "", rownames(rscu.grid)), metadata$file)]
#Remove outgroup
rscu.grid <- rscu.grid[match(metadata$name[metadata$ingroup == 1], rownames(rscu.grid)),]
#Remove codons
rscu.grid <- rscu.grid[, !(colnames(rscu.grid) %in% c("tga", "tgg", "tag", "taa", "atg"))]

#Normalise data
rscu.grid <- scale(rscu.grid)
#Make distance matrix
rscu.dist <- dist(rscu.grid, method="euclidean")
#Do hierarchical clustering
rscu.hclust <- hclust(rscu.dist, method="average")

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
                       guide=guide_colourbar(title="Codon usage bias\n(RSCU)",
                                             title.position="top",
                                             direction="horizontal",)) +
  theme(legend.position=c(0.2, 0.8))

## COMBINE ##

#Write to file
tiff(file=paste0("Fig5-", Sys.Date(), ".tiff"),
     height=4, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(plot_grid(gg.pca,
                    gg.pcatree +
                      annotation_custom(ggplotGrob(gg.codongenus), ymin=55, ymax=65, xmin=-25, xmax=-45),
                    rel_widths=c(1.5, 1), align="h", axis="tb",
                    labels=c("A", ""), label_size=10),
          gg.rscu,
          labels="AUTO", label_size=10, ncol=2, rel_widths=c(0.85, 1))
dev.off()