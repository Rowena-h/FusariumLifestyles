library(plyr)
library(dplyr)
library(stringr)
library(ape)
library(ggtree)
library(ggthemes)
library(ggalluvial)
library(ggfittext)
library(UpSetR)
library(ggplotify)
library(ggpubr)
library(vegan)
library(ggnewscale)
library(tidyr)
library(ggupset)
library(colorspace)
library(jsonlite)
library(scales)
library(patchwork)
library(seqmagick)
library(ggmsa)
library(multcompView)
library(ggstance)
library(Biostrings)
library(dendextend)
library(matrixStats)
library(seqinr)
library(grid)
library(pBrackets)
library(reshape2)
library(MCMCtreeR)
library(deeptime)
library(phytools)
library(ggrepel)
library(eulerr)

#Colour palette
show_col(colorblind_pal()(8))

#Read in orthogroup data
load("effector_prediction/effector-matrices-2021-05-25.RData")

#Read in sample metadata
metadata <- read.csv("metadata.csv")

#Make dataframe of lifestyle colours
col.df <- data.frame(lifestyle=c("endophyte", "coral associated", "human pathogen", "insect associated", "plant associated", "plant pathogen", "saprotroph", "mycoparasite"),
                     colour=c("#009E73", "#FFE983", "#000000", "#56B4E9", "#9AE324", "dimgrey", "#0072B2", "#D55E00"))


##EFFECTOR PREDICTION SANKEY

protein.lengths <- list()

for (i in metadata$file2[metadata$ingroup != "outgroup"]) {
  
  #print(i)
  
  #Read in CSEPfilter log
  CSEPfilter <- read.csv(paste0("effector_prediction/CSEPfilter_", i, ".faa.log"), sep=":", row.names=NULL, header=FALSE)
  CSEPfilter <- CSEPfilter[!is.na(CSEPfilter$V2),]
  CSEPfilter$V1 <- word(CSEPfilter$V1, 1)
  CSEPfilter$V1[max(grep("Phobius", CSEPfilter$V1))] <- "Phobius2"
  rownames(CSEPfilter) <-CSEPfilter$V1
  CSEPfilter <- subset(CSEPfilter, select="V2")
  colnames(CSEPfilter) <- metadata$name[metadata$file2 == i]
  
  assign(paste0(i, ".CSEPfilter"), CSEPfilter)
  
  #Protein length stuff
  
  proteins <- read.fasta(paste0("orthology_inference/", i, ".faa"), seqtype="AA")
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
  
  effectors <- scan(paste0("effector_prediction/", i, ".faa_candidate_effectors"), character(), quote="")
  assign(paste0(i, ".effectors"), effectors)
  
}

CSEPfilter.df <- bind_cols(mget(paste0(metadata$file2[metadata$ingroup != "outgroup"], ".CSEPfilter")))
CSEPfilter.labels <- data.frame(programme=c("", "SignalP v5.0b", "TargetP v2.0", "Phobius v1.01", "TMHMM v2.0c", "Phobius v1.01", "ps_scan v1.86", "NucPred v1.1",  "PredGPI", "EffectorP v3.0"),
                              x=rev(seq(1.55, length(rownames(CSEPfilter.df)) + 1)),
                              mean=round(rowMeans(CSEPfilter.df)),
                              range=paste0(prettyNum(rowMins(as.matrix(CSEPfilter.df)), big.mark=","),
                                           " - ", 
                                           prettyNum(rowMaxs(as.matrix(CSEPfilter.df)), big.mark=",")),
                              sum=rowSums(CSEPfilter.df),
                              arrow.x=c(seq(2, 10), NA),
                              arrow.xend=c(seq(1.2, 9.2), NA),
                              face=c("bold", "plain", "plain", "bold", "plain", "plain", "plain", "plain", "plain", "bold"))
#CSEPfilter.df[1,] <- CSEPfilter.df[1,] / 5

CSEPfilter.df$programme <- rownames(CSEPfilter.df)
CSEPfilter.df <- melt(CSEPfilter.df)
CSEPfilter.df$programme <- factor(CSEPfilter.df$programme, levels=c("Total", "SignalP", "TargetP", "Phobius", "TMHMM", "Phobius2", "Prosite", "NucPred", "PredGPI", "EffectorP"))

gg.CSEPfilter <- ggplot(CSEPfilter.df,
                      aes(x=programme, y=value, stratum=variable, alluvium=variable, fill=variable)) +
  geom_stratum(colour=NA, size=0.2) +
  geom_flow(size=0.2) +
  geom_segment(data=CSEPfilter.labels,
               aes(x=arrow.x, xend=arrow.xend, y=-130000, yend=-130000),
               arrow=arrow(length=unit(0.1, "cm"), type="closed"),
               size=0.15,
               inherit.aes=FALSE) +
  geom_label(data=CSEPfilter.labels,
            aes(x=rev(c(1:length(rownames(CSEPfilter.labels)))), y=-130000, label=range),
            size=2.5,
            inherit.aes=FALSE) +
  geom_label(data=CSEPfilter.labels[CSEPfilter.labels$programme != "",],
             aes(x=x, y=-130000, label=programme),
             label.padding=unit(0.15, "lines"),
             label.size=NA,
             size=2,
             inherit.aes=FALSE) +
  annotate("text", x=10.5, y=-130000, label="Number of proteins", size=3) +
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
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, face=rev(CSEPfilter.labels$face)),
        axis.text.x=element_blank(),
        axis.title=element_blank())

#https://stackoverflow.com/questions/35633239/add-curly-braces-to-ggplot2-and-then-use-ggsave
bracketsGrob <- function(...){
  l <- list(...)
  e <- new.env()
  e$l <- l
  grid:::recordGrob(  {
    do.call(grid.brackets, l)
  }, e)
}

b1 <- bracketsGrob(0.025, 0.68, 0.025, 0.83, lwd=1.5, ticks=0.8, col="dimgrey", h=0.02)
b2 <- bracketsGrob(0.025, 0.38, 0.025, 0.53, lwd=1.5, ticks=0.8, col="dimgrey", h=0.02)

gg.CSEPfilter <- gg.CSEPfilter +
  annotation_custom(b1) +
  annotation_custom(b2)


lengths.df <- stack(protein.lengths)

effectors <- unlist(mget(paste0(metadata$file2[metadata$ingroup != "outgroup"], ".effectors")))

lengths.df$group[!is.na(match(rownames(lengths.df), effectors))] <- "effector"
lengths.df$group[is.na(lengths.df$group)] <- "other"

lengths.labels <- data.frame(table(lengths.df$group))
lengths.labels$mean <- t.test(lengths.df$values[lengths.df$group == "effector"], lengths.df$values[lengths.df$group == "other"])$estimate
lengths.labels$label <- paste0("n=", prettyNum(lengths.labels$Freq, big.mark=","))

gg.lengths <- ggplot(lengths.df, aes(x=group, y=values)) +
  geom_hline(yintercept=300, lty="dashed") +
  geom_violin(size=0.3) +
  geom_boxplot(width=0.2, size=0.3, outlier.size=0.3) +
  geom_text(data=lengths.labels, aes(x=Var1, y=mean, label=label),
            nudge_y=500, nudge_x=0.3, size=1.5, inherit.aes=FALSE) +
  scale_x_discrete(labels=c("CSEP", "Other"),
                   limits=c("effector", "other")) +
  scale_y_continuous(labels=comma, sec.axis=dup_axis(),
                     breaks=c(300, 5000, 10000),
                     #expand=c(0, 100),
                     minor_breaks=seq(0, round_any(max(lengths.df$values), 1000, f=ceiling), 1000)) +
  labs(y="Protein length") +
  theme_minimal() +
  theme(panel.grid.major.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=15),
        axis.title.y.left=element_text(size=15),
        axis.title.y.right=element_blank(),
        axis.text.y.left=element_blank(),
        axis.text.y.right=element_text(size=10),
        panel.border=element_rect(colour="black", fill=NA, size=1.5))

tiff(file=paste0("CSEP-prediction-plot-", Sys.Date(), ".tiff"), width=8, height=5, units="in", res=300, compression="lzw")
gg.CSEPfilter + 
  annotation_custom(ggplotGrob(gg.lengths), ymin=200000, ymax=900000, xmin=2, xmax=8)
dev.off()


#ENDOPHYTE DETERMINANTS

enriched.df <- read.csv("endophyte_determinants_test/retainedOGs.csv")

orthogroups.stats[match(enriched.df$ogs[enriched.df$enrichedInEndophytes == "True"], orthogroups.stats$orthogroup),]




##TANGLEGRAM - SUPPLEMENTARY FIG 1##

astral <- read.tree("phylogenomics/species_tree/astral/fus_astral_proteins_62T.tre")
astral$edge.length <- rep(1, length(astral$edge.length))
raxmlng <- read.tree("phylogenomics/species_tree/raxml-ng/fus_proteins_62T.raxml.support")
iqtree <- read.tree("phylogenomics/species_tree/iqtree/fus_proteins_62T_iqtree_genepart.contree")

#Proportion of supported branches (ultrafast bootstrap >= 95, standard bootstrap >= 70)
suppressWarnings(round(length(which(as.numeric(raxmlng$node.label) >= 70)) /
                   length(na.omit(as.numeric(raxmlng$node.label))) * 100))
suppressWarnings(round(length(which(as.numeric(iqtree$node.label) >= 95)) /
                   length(na.omit(as.numeric(iqtree$node.label))) * 100))


outgroup <- "Ilyonectria sp."

for (i in c("astral", "raxmlng", "iqtree")) {
  tree <- get(i)
  tree$tip.label <- metadata$name[match(tree$tip.label, metadata$tip)]
  tree <- root(tree, outgroup, resolve.root=TRUE, edgelabel=TRUE)
  tree <- suppressWarnings(chronos(tree))
  tree <- as.dendrogram(tree)
  assign(paste0(i, ".dend"), tree)
}

#Untangle dendrograms
dend <- untangle_labels(raxmlng.dend, astral.dend, method="step2side")

#line.col <- rep("grey", length(labels(iqtree.dend)))
#line.col[c(59:62)] <- "black"

#tiff(file=paste0("tanglegram-", Sys.Date(), ".tiff"), height=8, width=8, units="in", res=300, compression="lzw")

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

#dev.off()


## MAIN TREE PLOT ##

for (i in c("iqtree", "raxmlng")) {
  
  #Root ML tree
  tree <- root(get(i), "Ilysp1_GeneCatalog_proteins_20121116", resolve.root=TRUE, edgelabel=TRUE)
  #Edit tip names
  tree$tip.label <- metadata$name[match(tree$tip.label, metadata$tip)]
  
  if (i == "iqtree") {
    tree$node.label[which(suppressWarnings(as.numeric(tree$node.label) < 95))] <- "∗"
    tree$node.label[tree$node.label != "∗"] <- ""
  } else if (i == "raxmlng") {
    tree$node.label[which(suppressWarnings(as.numeric(tree$node.label) < 70))] <- "†"
    tree$node.label[tree$node.label != "†"] <- ""
  }
  
  #Plot ML tree
  gg <- ggtree(tree, branch.length="none") %<+% metadata
  #Capture ML data structure
  gg.tree.data <- gg[["data"]] %>%
    arrange(y)
  
  assign(paste0(i, ".tree"), tree)
  assign(paste0("gg.", i), gg)
  assign(paste0("gg.tree.data.", i), gg.tree.data)
  
}

ml.nodes <- matchNodes(raxmlng.tree, iqtree.tree)

species.tree <- iqtree.tree
species.tree$node.label <- paste0(species.tree$node.label, gg.tree.data.raxmlng$label[match(ml.nodes[,1], gg.tree.data.raxmlng$node)][order(ml.nodes[,2])])

#Plot ML tree
gg.speciestree <- ggtree(species.tree, branch.length="none") %<+% metadata

#Capture ML data structure
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
  dated.tree$apePhy$node.label <- gg.tree.data.speciestree$label[match(nodes[,1], gg.tree.data.speciestree$node)][order(nodes[,2])]
  
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
  tiplabel.face.datedtree[which(metadata$new[match(dated.tree$apePhy$tip.label, metadata$name)] == "Y")] <- "bold.italic"
  
  #Make dataframe for plotting divergence time confidence intervals
  confidence.intervals <- as.data.frame(dated.tree$nodeAges)
  confidence.intervals <- confidence.intervals * 100
  confidence.intervals$ymin <- gg.tree.data.dated$y[match(rownames(confidence.intervals), gg.tree.data.dated$node)] -0.2
  confidence.intervals$ymax <- gg.tree.data.dated$y[match(rownames(confidence.intervals), gg.tree.data.dated$node)] +0.2
  
  #Add confidence intervals and time scale to tree plot
  gg.datedtree <- gg.datedtree +
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
    theme(plot.margin=unit(c(1, 0, 0, 2), "cm"),
          axis.text.x.bottom=element_text(size=5),
          axis.title.x.bottom=element_text(size=5))
  
  #Reverse the x axis
  gg.datedtree <- revts(gg.datedtree)

  #Make dataframe of species complex nodes
  sc.df.dated <- data.frame(sc=names(table(metadata$speciescomplex.abb[match(dated.tree$apePhy$tip.label, metadata$name)]))[table(metadata$speciescomplex.abb[match(dated.tree$apePhy$tip.label, metadata$name)]) > 1],
                            node=NA)
  
  #Get nodes for each species complex
  for (i in 1:length(sc.df.dated$sc)) {
    sc.df.dated$node[i] <- getMRCA(dated.tree$apePhy, metadata$name[metadata$speciescomplex.abb == sc.df.dated$sc[i]])
  }
  
  #Make dataframe for species complexes with only one taxon represented
  sc.df.dated2 <- data.frame(sc=unique(metadata$speciescomplex.abb)[is.na(match(unique(metadata$speciescomplex.abb), sc.df.dated$sc))],
                             node=NA)
  
  sc.df.dated2$node <- gg.tree.data.dated$node[match(sc.df.dated2$sc, gg.tree.data.dated$speciescomplex.abb)]
  
  #Combine dataframes
  sc.df.dated <- rbind(sc.df.dated, sc.df.dated2)
  sc.df.dated <- sc.df.dated[match(na.omit(unique(gg.tree.data.dated$speciescomplex.abb)), sc.df.dated$sc),]
  #Make alternated coding for highlights on tree
  sc.df.dated$box <- rep(c(0,1), length.out=length(sc.df.dated$sc))
  
  #Add species complex highlights and labels to tree plot
  for (i in 1:length(sc.df.dated$sc)) {
    
    if (sc.df.dated$box[i] == 1) {
      
      gg.datedtree <- gg.datedtree +
        geom_highlight(node=sc.df.dated$node[sc.df.dated$sc == sc.df.dated$sc[i]],
                       fill="#000000", alpha=0.075, extend=200)
      
    } else {
      
      gg.datedtree <- gg.datedtree +
        geom_highlight(node=sc.df.dated$node[sc.df.dated$sc == sc.df.dated$sc[i]],
                       fill="#000000", alpha=0.04, extend=200)
      
    }
    
    gg.datedtree <- gg.datedtree +
      geom_cladelabel(node=sc.df.dated$node[sc.df.dated$sc == sc.df.dated$sc[i]],
                      label=sc.df.dated$sc[i], offset.text=2, offset=60, align=TRUE, barsize=0.2, fontsize=1.3)
    
  }
  
  #Add final annotations to tree plot
  gg.datedtree <- gg.datedtree +
    geom_vline(xintercept=epochs$max_age[which((epochs$max_age) < max(confidence.intervals$`95%_upper`))] * -1,
               linetype="dashed", colour="grey", size=0.2) +
    geom_nodelab(size=1.4, hjust=1.5, vjust=-0.5) +
    geom_tiplab(aes(label=tiplab), fontface=tiplabel.face.datedtree, size=1.5, offset=1) +
    geom_tree(size=0.25) +
    geom_tippoint(aes(colour=lifestyle), size=0.7, show.legend=FALSE) +
    geom_nodepoint(aes(subset=(node %in% c(64, 66))),
                   size=1, colour="black") +
    geom_text_repel(data=gg.tree.data.dated[gg.tree.data.dated$node %in% c(64, 66),],
                    aes(x=x-max(gg.tree.data.dated$x), y=y,
                        label=c("O'Donnell et al. 2020\ngeneric concept", "Lombard et al. 2015\ngeneric concept")),
                    hjust=1,
                    size=1.1,
                    segment.size=0.25,
                    min.segment.length=unit(0, 'lines'),
                    xlim=c(-Inf, Inf),
                    nudge_x=-3,
                    nudge_y=2.5) +
    scale_colour_manual(values=col.df$colour[match(sort(unique(metadata$lifestyle[match(dated.tree$apePhy$tip.label,
                                                                                        metadata$name)])),
                                                   col.df$lifestyle)]) +
    scale_linetype_manual(values=c("solid", "11"))
  
  assign(paste0("gg.datedtree.", clock), gg.datedtree)
  assign(paste0("gg.tree.data.dated.", clock), gg.tree.data.dated)
  
}

#Write to file
tiff(file=paste0("dated-trees-", Sys.Date(), ".tiff"),
     height=8, width=6, unit="in", res=600, compression="lzw")
plot_grid(gg.datedtree.correlated, gg.datedtree.independent, ncol=1, labels="AUTO", align="v", axis="lr")
dev.off()


##CORE ACCESSORY SPECIFIC BARGRAPH

orthogroups.copies <- get("orthogroups.copies.ingroup2")
orthogroups.stats <- get("orthogroups.stats.ingroup2")

secretome.df <- data.frame(taxon=rep(colnames(orthogroups.copies), each=3), 
                           secretome=rep(c("specific", "accessory", "core"), length(colnames(orthogroups.copies))),
                           orthogroups=NA,
                           effectors=NA)

for (i in unique(secretome.df$taxon)) {
  for (j in unique(secretome.df$secretome)) {
    
    secretome.df$orthogroups[intersect(which(secretome.df$taxon == i), which(secretome.df$secretome == j))] <- table(orthogroups.stats$secretome[match(rownames(orthogroups.copies[orthogroups.copies[, i] > 0,]), orthogroups.stats$orthogroup)])[j]
    
    secretome.df$effectors[intersect(which(secretome.df$taxon == i), which(secretome.df$secretome == j))] <- table(orthogroups.stats$secretome[match(rownames(effector.count[effector.count[, i] > 0,]), orthogroups.stats$orthogroup)])[j]
    
  }
}

secretome.df$taxon <- metadata$name[match(secretome.df$taxon, metadata$file2)]
secretome.df$secretome <- factor(secretome.df$secretome, levels=c("specific", "accessory", "core"))
secretome.df$taxon <- factor(secretome.df$taxon, levels=rev(tip.order.datedtree))

secretome.df$y <- as.numeric(factor(secretome.df$taxon))


#secretome.df <- melt(secretome.df)

gg.secretome.orthogroups <- ggplot(secretome.df, aes(y=taxon, x=orthogroups, fill=secretome)) +
  geom_bar(stat="identity", size=0.5, width=0.6) +
  labs(x="Orthogroups") +
  scale_y_discrete(limits=rev(tip.order.datedtree),
                   expand=expansion(0, 0)) +
  coord_cartesian(clip="off",
                  ylim=c(0, 62)) +
  scale_x_continuous(expand=c(0, 0),
                     position="top",
                     labels=scales::comma) +
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
        panel.grid.major.y=element_blank())

gg.secretome.effectors <- ggplot(secretome.df, aes(y=taxon, x=effectors, fill=secretome)) +
  geom_bar(stat="identity", size=0.5, width=0.6) +
  labs(x="CSEPs") +
  scale_y_discrete(limits=rev(tip.order.datedtree),
                   expand=expansion(0, 0)) +
  coord_cartesian(clip="off",
                  ylim=c(0, 62)) +
  guides(fill=guide_legend(nrow=1)) +
  scale_x_continuous(expand=c(0, 0),
                     position="top") +
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
        legend.text=element_text(size=5),
        legend.key.size=unit(0.2, "cm"),
        legend.spacing.x=unit(0.1, "cm"),
        legend.title=element_blank())


lifestyle.grid <- metadata[c(14:21)]
rownames(lifestyle.grid) <- metadata$name

##Fix colours

col.df$lifestyle <- factor(col.df$lifestyle, levels=c("plant associated", "endophyte", "plant pathogen", "saprotroph", "insect associated", "mycoparasite", "human pathogen", "coral associated"))

gg.lifestyles.grid <- gheatmap(gg.datedtree.correlated,
                               lifestyle.grid,
                               offset=73,
                               width=0.15, 
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
  annotate("text", x=80, y=65, label="All reported\nlifestyles", size=2) +
  theme(plot.margin=unit(c(0, 1, 0, 0), "cm"),
        legend.position="top",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(5, "pt"),
        legend.text=element_text(size=4, margin=margin(0, 5, 0, 0)),
        legend.title=element_text(face="bold", size=5, margin=margin(0, 5, 0, 0)))


tiff(file=paste0("lifestyles-tree-", Sys.Date(), ".tiff"),
     height=4, width=6.75, unit="in", res=600, compression="lzw")
plot_grid(gg.lifestyles.grid, gg.secretome.orthogroups, gg.secretome.effectors,
          nrow=1, rel_widths=c(5, 1, 1), align="h", axis="bt")
dev.off()


##SPECIFIC EFFECTORS

orthogroups.count.specific <- orthogroups.copies.ingroup2[!is.na(match(rownames(orthogroups.copies.ingroup2), orthogroups.stats.ingroup2$orthogroup[intersect(which(is.na(orthogroups.stats.ingroup2$effector)),                                                            which(orthogroups.stats.ingroup2$secretome == "specific"))])),]

effector.count.specific <- effector.count.ingroup2[!is.na(match(rownames(effector.count.ingroup2), orthogroups.stats.ingroup2$orthogroup[orthogroups.stats.ingroup2$secretome == "specific"])),]

specific.df <- data.frame(species=colnames(effector.count.specific),
                          effectors=NA,
                          other=NA)

for (i in specific.df$species) {
  specific.df$effectors[specific.df$species == i] <- length(which(effector.count.specific[,i] > 0))
  specific.df$other[specific.df$species == i] <- length(which(orthogroups.count.specific[,i] > 0))
  
}

specific.df$lifestyle <- metadata$lifestyle[match(specific.df$species, metadata$file2)]
specific.df$lifestyle <- sub("-", " ", specific.df$lifestyle)

#Tukey significance testing
tukey <- TukeyHSD(aov(lm(effectors ~ lifestyle, data=specific.df)))
#Make dataframe for ggplot with tukey groups
tukey.df <- data.frame(multcompLetters(tukey[["lifestyle"]][,4])["Letters"])
tukey.df <- data.frame(Treatment=rownames(tukey.df), Letters=tukey.df$Letters)

specific.df <- melt(specific.df)

ggplot(specific.df, aes(x=lifestyle, y=value, fill=lifestyle)) +
  facet_wrap(. ~ variable, scales="free") +
  geom_boxplot() +
  #geom_text(data=tukey.df,
  #          aes(x=Treatment, y=Inf, label=Letters),
  #          family="mono",
  #          hjust=0.5,
  #          size=3,
  #          inherit.aes=FALSE) +
  scale_fill_manual(values=col.df$colour[match(sort(unique(specific.df$lifestyle)), col.df$lifestyle)],
                    labels=str_to_sentence(col.df$lifestyle[match(sort(unique(specific.df$lifestyle)),
                                                                  col.df$lifestyle)]),
                    name=NULL) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())



orthogroups.count.accessory <- as.data.frame(lapply(orthogroups.copies.ingroup2[!is.na(match(rownames(orthogroups.copies.ingroup2), orthogroups.stats.ingroup2$orthogroup[orthogroups.stats.ingroup2$secretome == "accessory"])),], as.logical))

accessory.df <- as.data.frame(table(rowSums(orthogroups.count.accessory)))
accessory.df$Var1 <- factor(accessory.df$Var1)

ggplot(accessory.df, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity")


##LIFESTYLE TEST RESULTS

#For effectors and orthogroups...
for (i in c("effectors", "orthogroups")){
  
  print(i)
  #Read in test results
  phy.pca.result <- read.csv(paste0("lifestyle_test/", i, "/metadata.csv"))
  lifestyle.data <- read.csv(paste0("lifestyle_test/", i, "/data.csv"), row.names="genome")
  #Make distance matrix
  dist <- vegdist(lifestyle.data, method="jaccard")
  #Do permanova
  permanova <- adonis2(formula=dist ~ PC1 + PC2 + lifestyle, data=phy.pca.result, permutations=9999)
  
  print(paste0("Phylogeny: ", round(sum(permanova$R2[1:2]) * 100), "%"))
  print(paste0("Lifestyle: ", round(sum(permanova$R2[3]) * 100), "%"))
  
  pca.result <- read.csv(paste0("lifestyle_test/", i, "/pca.csv"))
  pca.result$taxon <- metadata$name[match(pca.result$genome, metadata$tip)]
  
  gg.pca <- ggplot(pca.result, aes(x=PC1, y=PC2, colour=lifestyle)) +
    geom_point() +
    geom_text_repel(aes(label=sub("Fusarium", "F.", taxon)), size=2, max.overlaps=20) +
    scale_colour_manual(values=col.df$colour[match(sort(unique(metadata$lifestyle)), col.df$lifestyle)],
                        labels=sort(unique(metadata$lifestyle)),
                        guide=FALSE) +
    theme_void() +
    theme(panel.border=element_rect(colour="black", size=2, fill=NA))
  
  assign(paste0("permanova.", i), permanova)
  assign(paste0("gg.pca.", i), gg.pca)
  
}


pw.lifestyle.test <- rbind(data.frame(melt(as.matrix(read.csv("lifestyle_test/effectors/pairwiseComparisons.csv",
                                                              row.names=1)), na.rm=TRUE), data="CSEPs"),
                           data.frame(melt(as.matrix(read.csv("lifestyle_test/orthogroups/pairwiseComparisons.csv",
                                                              row.names=1)), na.rm=TRUE), data="Orthogroups"))

pw.lifestyle.test$value <- round(pw.lifestyle.test$value, digits=3)
pw.lifestyle.test$Var2 <- sub("p.value.", "", pw.lifestyle.test$Var2)

#Plot grid
gg.pwperm <- ggplot(pw.lifestyle.test, aes(Var2, Var1, fill=value>0.05)) +
  facet_grid(. ~ data, switch="x") +
  geom_tile(color="grey", size=2, show.legend=FALSE) +
  geom_text(aes(label=value), size=3) +
  scale_fill_manual(values=c("white", "darkgrey")) +
  scale_x_discrete(position="top") + 
  theme_minimal() + 
  theme(axis.text.x=element_text(colour=c("#009E73","#56B4E9", "#D55E00", "#9AE324", "dimgrey"), size=4),
        axis.text.y=element_text(colour=c("#56B4E9", "#D55E00", "#9AE324", "dimgrey", "#0072B2"),
                                 angle=90, hjust=0.5, size=4),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text=element_text(face="bold", size=9),
        panel.grid=element_blank(),
        plot.margin=unit(c(0, 10, 0, 10), "pt")) +
  coord_fixed()

pw.lifestyle.test$Var1 <- as.character(pw.lifestyle.test$Var1)

pw.lifestyle.test$value[pw.lifestyle.test$value <= 0.05] <- 0

network.CSEP <- network(pw.lifestyle.test[pw.lifestyle.test$data == "CSEPs",], matrix.type="edgelist")

network.CSEP %v% "lifestyle" = network.vertex.names(network.CSEP)

cols <- col.df$colour[match(sort(unique(metadata$lifestyle)), col.df$lifestyle)]
names(cols) <- sub(" ", "", sort(unique(metadata$lifestyle)))

ggnet2(network.CSEP,
       color="lifestyle",
       palette=cols,
       label=TRUE)

tiff(file=paste0("pvalue-grid-", Sys.Date(), ".tiff"),
     width=6, height=3, units="in", res=300, compression="lzw")
gg.pwperm
dev.off()



##SELECTION

#Vector of core, single-copy orthogroups
core.SC.orthogroups <- Reduce(intersect,
                              list(orthogroups.stats.ingroup0$orthogroup[which(
                                orthogroups.stats.ingroup0$copy_number == "single")],
                                orthogroups.stats.ingroup0$orthogroup[which(
                                  orthogroups.stats.ingroup0$secretome == "core")]))

#Make dataframe for selection test results
selection.df <- data.frame(orthogroup=core.SC.orthogroups,
                           effector="N", busted="N", absrel="N", busted.success="Y", absrel.success="Y")
#Orthogroups predicted to be effectors
selection.df$effector[match(core.SC.mixed, selection.df$orthogroup)] <- "Y"

#Make empty list to capture aBSREL p-values for each lineage
absrel.p <- list()

#Create bar to show progress
progress.bar <- txtProgressBar(1, length(selection.df$orthogroup), initial=0, char="=", style=3)

#For each orthogroup...
for (i in selection.df$orthogroup) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, which(selection.df$orthogroup == i))
  
  #Try to read in BUSTED results and add whether the p-value is significant to the results dataframe
  busted.results <- tryCatch(fromJSON(paste0("selection/hyphy/busted/", i, "_BUSTED.json")), error=function(e) NULL)
  
  if (!is.null(busted.results)) {
    
    if (busted.results$`test results`$`p-value` < 0.05)
      
      selection.df$busted[match(i, selection.df$orthogroup)] <- "Y"
    
  } else {
    selection.df$busted.success[match(i, selection.df$orthogroup)] <- "N"
  }
  
  #Try to read in aBSREL results and add p-value results to list
  absrel.results <- tryCatch(fromJSON(paste0("selection/hyphy/absrel/", i, "_aBSREL.json")), error=function(e) NULL)
  
  if (!is.null(absrel.results)) {
    
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

#Make matrix for euler
selection.mat.labels <- c("Total"=length(union(which(selection.df$busted.success == "Y"),
                                               which(selection.df$absrel.success == "Y"))),
                          "BUSTED"=0,
                          "aBSREL"=0,
                          "CSEPs"=0,
                          "BUSTED&aBSREL"=0,
                          "CSEPs&BUSTED"=0,
                          "CSEPs&aBSREL"=0,
                          "Total&BUSTED"=length(Reduce(intersect,
                                                       list(which(selection.df$effector == "N"),
                                                            which(selection.df$busted == "Y"),
                                                            which(selection.df$absrel == "N")))),
                          "Total&aBSREL"=length(Reduce(intersect,
                                                       list(which(selection.df$effector == "N"),
                                                            which(selection.df$busted == "N"),
                                                            which(selection.df$absrel == "Y")))),
                          "Total&CSEPs"=length(Reduce(intersect,
                                                      list(which(selection.df$effector == "Y"),
                                                           which(selection.df$busted == "N"),
                                                           which(selection.df$absrel == "N")))),
                          "Total&BUSTED&aBSREL"=length(Reduce(intersect,
                                                              list(which(selection.df$effector == "N"),
                                                                   which(selection.df$busted == "Y"),
                                                                   which(selection.df$absrel == "Y")))),
                          "Total&BUSTED&CSEPs"=length(Reduce(intersect,
                                                             list(which(selection.df$effector == "Y"),
                                                                  which(selection.df$busted == "Y"),
                                                                  which(selection.df$absrel == "N")))),
                          "Total&aBSREL&CSEPs"=length(Reduce(intersect,
                                                             list(which(selection.df$effector == "Y"),
                                                                  which(selection.df$busted == "N"),
                                                                  which(selection.df$absrel == "Y")))),
                          "Total&BUSTED&aBSREL&CSEPs"=length(Reduce(intersect,
                                                                    list(which(selection.df$effector == "Y"),
                                                                         which(selection.df$busted == "Y"),
                                                                         which(selection.df$absrel == "Y")))))

print(paste0(selection.mat.labels[["Total&BUSTED&aBSREL&CSEPs"]] + selection.mat.labels[["Total&BUSTED&aBSREL"]],
      ", ",
      (selection.mat.labels[["Total&BUSTED&aBSREL&CSEPs"]] + selection.mat.labels[["Total&BUSTED&aBSREL"]]) /
        selection.mat.labels[["Total"]],
      "%"))

selection.mat <- selection.mat.labels
selection.mat[["Total&BUSTED&CSEPs"]] <- 9
#selection.mat[["Total&aBSREL&CSEPs"]] <- 8
#selection.mat.labels[["Total&BUSTED&aBSREL&CSEPs"]] <- 11
selection.mat[["Total&CSEPs"]] <- 11

#Create and plot euler diagram
set.seed(1)
selection.euler <- euler(selection.mat, shape="ellipse")
eulerr_options(padding=unit(0.5, "pt"))
euler <- plot(selection.euler,
              fills=list(fill=c("white", "dimgrey", "dimgrey", "white", rep("", 6),
                                "#E69F00", "grey", "grey", rep("", 1), "#E69F00"),
                         alpha=0.3),
              labels=list(labels=c("Orthogroups",
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
                                       paste0("Positively\nselected\nconsensus\n",
                                              selection.mat.labels[["Total&BUSTED&aBSREL"]]), 
                                       selection.mat.labels[["Total&BUSTED&CSEPs"]],
                                       selection.mat.labels[["Total&aBSREL&CSEPs"]], "",
                                       selection.mat.labels[["Total&BUSTED&aBSREL&CSEPs"]]),
                              col=c(rep("dimgrey", 4), "black", rep("dimgrey", 2), "black"),
                              cex=0.3))

#Convert to grob
euler.grob <- as.grob(euler)

#Read in dataframe matching aBSREL tree nodes to IQ-TREE species tree nodes
absrel.df <- read.csv("selection/hyphy/absrel/absrel_nodes.csv")
#Add tip labels
absrel.df <- rbind(absrel.df, data.frame(node=1:length(iqtree$tip.label), absrel=iqtree$tip.label))

##Add number of significant aBSREL tests per taxon that are also supported by busted
absrel.p.count <- lengths(lapply(absrel.p.sig, function(x)
  x[which(x %in% selection.df$orthogroup[which(selection.df$busted == "Y")])]))
absrel.df$num <- absrel.p.count[match(absrel.df$absrel, names(absrel.p.count))]

#Add names of positively selected effectors
absrel.p.names <- unlist(lapply(absrel.p.sig, function(x)
  paste(x[intersect(which(x %in% core.SC.mixed),
                    which(x %in% selection.df$orthogroup[which(selection.df$busted == "Y")]))],
          collapse=" ")))
absrel.p.names <- absrel.p.names[absrel.p.names != ""]
absrel.p.names <- gsub("OG000", "", absrel.p.names)

absrel.df$effectors <- NA

absrel.df$effectors <- absrel.p.names[match(absrel.df$absrel, names(absrel.p.names))]
absrel.df$effectors[absrel.df$effectors == ""] <- NA
absrel.df$effectors <- gsub('(?=(?:.{10})+$)', "\n", absrel.df$effectors, perl = TRUE)

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

#Plot fully annotated selection tree
gg.selection2 <- gg.selection %<+% absrel.df +
  geom_tree(size=1, aes(colour=num)) +
  xlim(0, 25) +
  scale_colour_gradient2(trans="pseudo_log",
                         low="black", mid="#F0E442", high="#CC79A7", midpoint=1,
                         breaks=c(0, 10, 100),
                         limits=c(0, 100),
                         labels=c(0 , 10, 100),
                         na.value="dimgrey",
                         guide=guide_colourbar(title="Orthogroups showing\npositive selection",
                                               title.position="top",
                                               direction="horizontal")) +
  new_scale_colour() +
  geom_tippoint(aes(colour=lifestyle), size=0.7) +
  scale_colour_manual(values=col.df$colour[na.omit(match(sort(unique(metadata$lifestyle[match(iqtree.tree$tip.label, metadata$name)])), col.df$lifestyle))],
                      labels=str_to_sentence(col.df$lifestyle[na.omit(match(sort(unique(metadata$lifestyle[match(iqtree.tree$tip.label, metadata$name)])), col.df$lifestyle))]),
                      na.translate=FALSE,
                      guide=guide_legend(title="Lifestyle",
                                         ncol=2,
                                         title.hjust=0)) +
  new_scale_colour() +
  geom_label2(aes(x=branch, label=effectors, fill=num),
              size=1, label.size=0, label.padding=unit(1.5, "pt"), fontface="plain") +
  scale_fill_gradient2(trans="pseudo_log",
                         low="black", mid="#F0E442", high="#CC79A7", midpoint=1,
                         breaks=c(0, 10, 100),
                         limits=c(0, 100),
                         labels=c(0 , 10, 100),
                         na.value="dimgrey",
                         guide=guide_colourbar(title="Orthogroups showing\npositive selection",
                                               title.position="top",
                                               direction="horizontal")) +
  geom_tiplab(aes(label=tiplab), fontface=tiplabel.face.speciestree, size=1.5, offset=0.1, colour="black") +
  new_scale_fill() +
  coord_cartesian(clip="off") +
  annotation_custom(euler.grob, xmin=-3, xmax=4, ymin=25, ymax=55) +
  theme(legend.box="horizontal",
        legend.title=element_text(size=5, face="bold"),
        legend.text=element_text(size=5),
        legend.position=c(0.15, 0.92),
        legend.key.size=unit(6, "pt"),
        legend.spacing=unit(0, "pt"),
        plot.margin=margin(0, 5.5, 0, 50, unit="pt"))


#Write to file
tiff(file=paste0("selectiontree-", Sys.Date(), ".tiff"),
     height=4, width=6.75, unit="in", res=600, compression="lzw")
gg.selection2
dev.off()


##MEME

msa.labels <- sub("[[:space:]]+\\(.*", "", metadata$tiplab)

meme.orthos <- list()

for (i in core.SC.mixed[which(core.SC.mixed %in% absrel.p.names)]) {
  
  meme.results <- fromJSON(paste0("selection/hyphy/meme/", i, "_MEME.json"))
  selection.sites <- which(meme.results$MLE$content$`0`[,7] <= 0.05)
  
  gene.tree <- read.tree(paste0("selection/trees/RAxML_bipartitions.", i, "_rooted"))
  gene.tree$tip.label <- sub("\\{FOREGROUND\\}", "", gene.tree$tip.label)
  #gene.tree$tip.label <- metadata$tip[match(gene.tree$tip.label, metadata$file)]
  
  #meme.orthos[i] <- i
  
  alignment <- fa_read(paste0("selection/alignments/codon/", i, "_aln_nuc.translated"))
  
  names(alignment) <- metadata$name[match(names(alignment), metadata$file2)]
  alignment <- alignment[match(tip.order.speciestree, names(alignment)),]
  names(alignment) <- sub("Fusarium", "F.", names(alignment))
  
  lines <- data.frame(x=0, xend=unique(width(alignment))+0.5, y=0:length(names(alignment))+0.5, yend=0:length(names(alignment))+0.5)
  rect <- data.frame(xmin=0, xmax=unique(width(alignment))+0.5, ymin=0.5, ymax=length(names(alignment))+0.5)
  sites <- data.frame(site=selection.sites)
  
  gg.msa <- ggmsa(alignment, font=NULL, border=FALSE, color="LETTER", posHighligthed=selection.sites) +
    geom_segment(data=lines, aes(x=x, xend=xend, y=y, yend=yend), colour="snow2", size=0.2) +
    geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, colour="dimgrey", size=0.2) +
    #geom_text_repel(data=sites, aes(label=site, x=site, y=29), 
    #                colour="dimgrey", size=2, angle=90, 
    #                #segment.color=NA,
    #                force=0.1) +
    ggtitle(sub("OG000", "", i)) +
    scale_x_continuous(breaks=selection.sites,
                       expand=c(0, 0),
                       #limits=c(0, unique(width(alignment))),
                       position="top") +
    coord_cartesian(clip="off", ylim=c(0, 61)) +
    theme(axis.text.x=element_text(size=3, hjust=0.5, vjust=0),
          axis.text.y=element_text(size=3),
          plot.title.position="plot",
          plot.title=element_text(size=7, face="bold", vjust=-5, hjust=0.3),
          plot.margin=margin(1, 1, 5, 1))  
  
  
  
  assign(paste0("alignment.", i), alignment)
  assign(paste0("gg.msa.", i), gg.msa)
  assign(paste0("meme.results.", i), meme.results)
  assign(paste0("selection.sites.", i), selection.sites)
  
}

#OG with no MEME sites
msa.to.plot <- sub("selection.sites.", "gg.msa.", names(mget(ls(pattern="selection.sites.")))[lengths(mget(ls(pattern="selection.sites."))) != 0])

gg.msa.gene <- do.call(ggarrange, c(mget(msa.to.plot), ncol=1))

tiff(file=paste0("OG0006864-alignment-", Sys.Date(), ".tiff"),
     height=3, width=2, unit="in", res=300, compression="lzw")

gg.msa.OG0006864

dev.off()


do.call(ggarrange, c(mget(paste0("gg.msa.gene.", rownames(effector.count.SC.SP[effector.count.SC.SP$secretome == "core",]))), ncol=1))




meme.df <- read.csv("meme_nodes.csv")
meme.results <- fromJSON("MEME_results/OG0005992_MEME.json")
meme.results.gene <- fromJSON("MEME_results/OG0005992_MEME_genetree.json")
selection.sites <- which(meme.results$MLE$content$`0`[,7] <= 0.05)
selection.sites.gene <- which(meme.results.gene$MLE$content$`0`[,7] <= 0.05)

meme.tree <- read.tree(text=paste0(meme.results$input$trees$`0`, ";"))
absrel.tree <- read.tree(text=paste0(absrel.results$input$trees$`0`, ";"))
meme.tree.gene <- read.tree(text=paste0(meme.results.gene$input$trees$`0`, ";"))




meme.ebf <- list()

for (i in selection.sites) {
  for (j in 1:length(names(meme.results[["branch attributes"]][["0"]]))) {
    meme.ebf[[names(meme.results[["branch attributes"]][["0"]][j])]] <- meme.results[["branch attributes"]][["0"]][[names(meme.results[["branch attributes"]][["0"]])[j]]][[paste0("EBF site ", i, " (partition 1)")]]
  }
  
  meme.df.tmp <- meme.df[c(1, 5)]
  meme.df.tmp <- rbind(meme.df.tmp, data.frame(node=1:length(species.tree$tip.label), species=species.tree$tip.label))
  meme.df.tmp$ebf <- unlist(meme.ebf)[match(meme.df.tmp$species, names(meme.ebf))]
  #meme.df.tmp <- meme.df.tmp[-which(is.na(meme.df.tmp$ebf)),]
  
  assign(paste0("meme.df.", i), meme.df.tmp)
}

species.tree$tip.label <- metadata$name[match(species.tree$tip.label, metadata$tip)]

for (i in selection.sites) {
  
  gg.meme <- ggtree(species.tree, branch.length="none", size=0.4) %<+% get(paste0("meme.df.", i)) +
    #geom_tiplab(fontface=tiplabel.face, size=1, colour=tiplabel.col) +
    xlim(c(0, 10)) +
    aes(colour=ebf/100 ^ log(0.5, base=0.01)) +
    scale_colour_gradient2(name="EBF",
                           low="black", mid="#F0E442", high="#CC79A7", midpoint=0.5,
                           breaks=c(0, 0.01, 1) ^ log(0.5, base=0.01),
                           limits=c(0, 1),
                           labels=c(0 , 1, ">=100"),
                           na.value="black",
                           oob=squish) +
    theme(legend.position="right",
          plot.margin=margin(0, 0, 0, 0))
  
  gg.meme <- gg.meme %<+% metadata +
    new_scale_color() +
    geom_tippoint(aes(colour=lifestyle.hyp1), size=1, show.legend=FALSE) +
    scale_colour_manual(values=c("#009E73", "#56B4E9", "dimgrey", "#0072B2"), 
                        labels=c("Endophyte", "Insect-mutualist", "Plant pathogen", "Saprotroph"), 
                        na.translate=FALSE)
  
  gg.meme <- flip(gg.meme, 30, 32)
  gg.meme <- flip(gg.meme, 39, 33)
  
  #assign(paste0("gg.meme.", i), gg.meme)
  
}

ggarrange(gg.meme.56, gg.meme.77, gg.meme.119, gg.meme.336, gg.meme.344, common.legend=TRUE)

meme.ebf.gene <- list()

for (i in selection.sites.gene) {
  for (j in 1:length(names(meme.results.gene[["branch attributes"]][["0"]]))) {
    meme.ebf.gene[[names(meme.results.gene[["branch attributes"]][["0"]][j])]] <- meme.results.gene[["branch attributes"]][["0"]][[names(meme.results.gene[["branch attributes"]][["0"]])[j]]][[paste0("EBF site ", i, " (partition 1)")]]
  }
  
  meme.df.tmp <- meme.df[c(1, 4)]
  meme.df.tmp <- rbind(meme.df.tmp, data.frame(node=1:length(gene.tree.5992$tip.label), gene=gene.tree.5992$tip.label))
  meme.df.tmp$ebf <- unlist(meme.ebf.gene)[match(meme.df.tmp$gene, names(meme.ebf.gene))]
  #meme.df.tmp <- meme.df.tmp[-which(is.na(meme.df.tmp$ebf)),]
  assign(paste0("meme.df.gene.", i), meme.df.tmp)
}

gene.tree.5992$tip.label <- metadata$name[match(gene.tree.5992$tip.label, metadata$tip)]

for (i in selection.sites.gene) {
  
  gg.meme <- ggtree(gene.tree.5992, branch.length="none", size=0.4) %<+% get(paste0("meme.df.gene.", i)) +
    #geom_tiplab(fontface=tiplabel.face, size=1, colour=tiplabel.col) +
    xlim(c(0, 10)) +
    aes(colour=ebf/100 ^ log(0.5, base=0.01)) +
    scale_colour_gradient2(name="EBF",
                           low="black", mid="#F0E442", high="#CC79A7", midpoint=0.5,
                           breaks=c(0, 0.01, 1) ^ log(0.5, base=0.01),
                           limits=c(0, 1),
                           labels=c(0 , 1, ">=100"),
                           na.value="black",
                           oob=squish) +
    theme(legend.position="right",
          plot.margin=margin(0, 0, 0, 0))
  
  gg.meme <- gg.meme %<+% metadata +
    new_scale_color() +
    geom_tippoint(aes(colour=lifestyle.hyp1), size=1, show.legend=FALSE) +
    scale_colour_manual(values=c("#009E73", "#56B4E9", "dimgrey", "#0072B2"), 
                        labels=c("Endophyte", "Insect-mutualist", "Plant pathogen", "Saprotroph"), 
                        na.translate=FALSE)
  
  gg.meme <- flip(gg.meme, 36, 42)
  gg.meme <- flip(gg.meme, 35, 43)
  
  assign(paste0("gg.meme.gene.", i), gg.meme)
  
}

ggarrange(gg.meme.gene.56, gg.meme.gene.77, gg.meme.gene.119, gg.meme.gene.336, gg.meme.gene.344, gg.meme.gene.200, common.legend=TRUE)


protein_sequences <- fa_read("MEME_results/OG0005992_aln_nuc.translated.fas")

names(protein_sequences) <- metadata$name[match(names(protein_sequences), metadata$file)]
protein_sequences <- protein_sequences[tip.order,]

gg.msa <- ggmsa(protein_sequences, font=NULL, border=FALSE, color="LETTER", posHighligthed=selection.sites.gene) +
  geom_hline(colour="white", yintercept=0:length(names(protein_sequences))+0.5, size=0.3) +
  scale_x_continuous(breaks=selection.sites.gene,
                     position="top") +
  coord_cartesian() +
  theme(panel.background=element_rect(fill="snow2", colour="white"),
        axis.text.x=element_text(size=4, angle=90, vjust=0.3, hjust=0),
        axis.text.y=element_text(size=4, face=label.face, colour=label.cols))

tiff(file=paste0("test-", Sys.Date(), ".tiff"), height=5, width=8, unit="in", res=300)

gg.msa / (((gg.meme.gene.56 | gg.meme.gene.77 | gg.meme.gene.119 | gg.meme.gene.200 | gg.meme.gene.336 | gg.meme.gene.344) / (gg.meme.56 | gg.meme.77 | gg.meme.119 | plot_spacer() | gg.meme.336 | gg.meme.344)) +
            plot_layout(guides="collect")) +
  plot_layout(heights=c(1, 2))

dev.off()


library(Biostrings)
x <- readAAStringSet("MEME_results/OG0005992_aln_nuc.translated")
names(x) <- metadata$name[match(names(x), metadata$file)]
d <- as.dist(stringDist(x, method = "hamming")/width(x)[1])
species.tree <- bionj(d)
data <- tidy_msa(x)

ggtree(species.tree, branch.length="none") +
  xlim_tree(20) +
  geom_tiplab() +
  geom_facet(geom=geom_msa, data=data, font=NULL, panel="msa", colour="LETTER", border=FALSE, posHighligthed=selection.sites.gene) +
  scale_x_continuous(breaks=selection.sites.gene) +
  #geom_seqlogo(panel="msa") +
  theme(legend.position="none",
        axis.text.x=element_text(colour="black"),
        strip.background=element_blank(),
        strip.text=element_blank(),
        plot.margin=margin(0, 0, 10, 0))




for (i in rownames(effector.count.SC.SP[effector.count.SC.SP$secretome == "core",])[-18]) {
  check <- ggarrange(
    ggtree(get(paste0("phy.selection.", i))) +
      xlim(0, 15) +
      geom_nodelab() +
      geom_tiplab(),
    ggtree(phy) +
      xlim(0, 15) +
      geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
      geom_tiplab()
  )
  plot(check)
}






for (i in selection.sites.gene) {
  gg.sites <- ggseqlogo(as.vector(protein_sequences), method = 'prob', col_scheme="clustalx") +
    scale_x_continuous(limits=c(i-2.5, i+2.5)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank())
  gg.sites <- gg.sites +
    annotate("rect", xmin=i-0.5, xmax=i+0.5, ymin=-0.01, ymax=1.01, alpha=0.1, fill="red")
  
  assign(paste0("gg.sites.", i), gg.sites)
  
  #plot(ggarrange(get(paste0("gg.sites.", i)),
  #          get(paste0("gg.meme.", i)), 
  #          get(paste0("gg.meme.gene.", i)),
  #          common.legend=TRUE, legend="bottom", ncol=1))
}

tiff(file=paste0("test-", Sys.Date(), ".tiff"), height=6, width=4, unit="in", res=300)
ggarrange(gg.sites.336, gg.meme.336, gg.meme.gene.336, common.legend=TRUE, legend="bottom", ncol=1)

dev.off()

for (i in ls(pattern="absrel.results.*")) {
  test <- get(i)
  test <- read.tree(text=paste0(test$input$trees[[1]], ";"))
  blah <- ggarrange(
    ggtree(test) +
      xlim(0, 15) +
      geom_tiplab() +
      geom_nodelab(),
    ggtree(phy) +
      xlim(0, 0.2) +
      geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
      geom_tiplab())
  plot(blah)
}


##STATS FOR TEXT##

#Number of effectors undergoing positive selection at some point in the genus
print(paste0(length(unique(unlist(lapply(absrel.p, function(x) paste(names(which(x <=0.05))))))), "/", length(rownames(effector.count.SC.SP[effector.count.SC.SP$secretome == "core",])[-18])))

