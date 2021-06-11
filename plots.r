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

#Colour palette
show_col(colorblind_pal()(8))

#Read in orthogroup data
load("effector_prediction/effector-matrices-2021-05-25.RData")
#Read in sample metadata
metadata <- read.csv("metadata.csv")


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




#TANGLEGRAM

astral <- read.tree("phylogenomics/species_tree/astral/fus_astral_proteins_62T.tre")
astral$edge.length <- rep(1, length(astral$edge.length))
iqtree <- read.tree("phylogenomics/species_tree/iqtree/gene_partitions/fus_proteins_62T_iqtree_genepart.contree")

outgroup <- "Ilyonectria sp."

for (i in c("astral", "iqtree")) {
  tree <- get(i)
  tree$tip.label <- metadata$name[match(tree$tip.label, metadata$tip)]
  tree <- root(tree, outgroup, resolve.root=TRUE, edgelabel=TRUE)
  tree <- chronos(tree)
  tree <- as.dendrogram(tree)
  assign(paste0(i, ".dend"), tree)
}

#Untangle dendrograms
dend <- untangle(astral.dend, iqtree.dend, method="step2side")

line.col <- rep("grey", length(labels(iqtree.dend)))
line.col[c(59:62)] <- "black"

tiff(file=paste0("tanglegram-", Sys.Date(), ".tiff"), height=8, width=8, units="in", res=300, compression="lzw")

tanglegram(dend,
           main_left="ASTRAL-III",
           main_right="RAxMLv8",
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
           #color_lines=line.col,
           rank_branches=TRUE)

dev.off()



##LIFESTYLE GRID

#Make dataframe of lifestyle colours
col.df <- data.frame(lifestyle=c("endophyte", "coral associated", "human pathogen", "insect associated", "plant associated", "plant pathogen", "saprotroph", "mycoparasite"),
                     colour=c("#009E73", "#FFE983", "#000000", "#56B4E9", "#9AE324", "dimgrey", "#0072B2", "#D55E00"))

dated.tree <- readMCMCtree("divergence_time_estimation/mcmctree/run1_independent/FigTree.tre", 
                           forceUltrametric=TRUE)

dated.tree$nodeAges <- dated.tree$nodeAges * 100
dated.tree$apePhy$edge.length <- dated.tree$apePhy$edge.length * 100

dated.tree$apePhy$tip.label <- metadata$name[match(dated.tree$apePhy$tip.label, metadata$short.tip)]

gg.datedtree <- ggtree(dated.tree$apePhy) %<+% metadata

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
#Fontface vector
tiplabel.face.datedtree <- rep("italic", length(dated.tree$apePhy$tip.label))
tiplabel.face.datedtree[which(metadata$new[match(dated.tree$apePhy$tip.label, metadata$name)] == "Y")] <- "bold.italic"

new.tiplabels.dated <- dated.tree$apePhy$tip.label

for (i in 1:length(new.tiplabels.dated)) {
  if (metadata$old.name[match(new.tiplabels.dated, metadata$name)][i] != "") {
    new.tiplabels.dated[i] <- paste0(new.tiplabels.dated[i], " (=", metadata$old.name[match(new.tiplabels.dated, metadata$name)][i], ")")
  }
}

new.tiplabels.dated <- sub("Fusarium", "F.", new.tiplabels.dated)


gg.datedtree <- gg.datedtree +
  coord_geo(clip="off", xlim=c(max(dated.tree$nodeAges[,3]) * -1, 70), ylim=c(-2,Ntip(dated.tree$apePhy)),
            dat=list("epochs", "periods"),
            pos=list("bottom", "bottom"),
            size=list(2,3),
            height=list(unit(0.5, "line"), unit(1, "line")),
            abbrv=list(TRUE, FALSE),
            lwd=0, alpha=0.5, neg=TRUE, expand=FALSE) +
  scale_x_continuous(breaks=seq(round(max(dated.tree$nodeAges[,3])) * -1, 0, 10), labels=rev(seq(0, round(max(dated.tree$nodeAges[,3])), 10)), name="Million years") +
  theme(plot.margin=unit(c(1, 0, 0, 2), "cm"),
        axis.text.x.bottom=element_text(),
        axis.title.x.bottom=element_text(size=8))

gg.datedtree <- revts(gg.datedtree)



#Make dataframe of species complex nodes
sc.df.dated <- data.frame(sc=names(table(metadata$speciescomplex.abb[match(dated.tree$apePhy$tip.label, metadata$name)]))[table(metadata$speciescomplex.abb[match(dated.tree$apePhy$tip.label, metadata$name)]) > 1],
                    node=NA)

#Get nodes for each species complex
for (i in 1:length(sc.df.dated$sc)) {
  sc.df.dated$node[i] <- getMRCA(dated.tree$apePhy, metadata$name[metadata$speciescomplex.abb == sc.df.dated$sc[i]])
}

sc.df.dated2 <- data.frame(sc=unique(metadata$speciescomplex.abb)[is.na(match(unique(metadata$speciescomplex.abb), sc.df.dated$sc))],
                     node=NA)

sc.df.dated2$node <- gg.tree.data.dated$node[match(sc.df.dated2$sc, gg.tree.data.dated$speciescomplex.abb)]

sc.df.dated <- rbind(sc.df.dated, sc.df.dated2)
sc.df.dated <- sc.df.dated[match(na.omit(unique(gg.tree.data.dated$speciescomplex.abb)), sc.df.dated$sc),]
sc.df.dated$box <- rep(c(0,1), length.out=length(sc.df.dated$sc))

for (i in 1:length(sc.df.dated$sc)) {
  
  if (sc.df.dated$box[i] == 1) {
    
    gg.datedtree <- gg.datedtree +
      geom_highlight(node=sc.df.dated$node[sc.df.dated$sc == sc.df.dated$sc[i]],
                     fill="#000000", alpha=0.05, extend=55)
    
  } else {
    
    gg.datedtree <- gg.datedtree +
      geom_highlight(node=sc.df.dated$node[sc.df.dated$sc == sc.df.dated$sc[i]],
                     fill="#000000", alpha=0.1, extend=55)
    
  }
  
  gg.datedtree <- gg.datedtree +
    geom_cladelabel(node=sc.df.dated$node[sc.df.dated$sc == sc.df.dated$sc[i]],
                    label=sc.df.dated$sc[i], offset.text=2, offset=55, align=TRUE, fontsize=3)
  
}

gg.datedtree <- gg.datedtree +
  geom_vline(xintercept=periods$max_age[which((periods$max_age) < max(dated.tree$nodeAges[,3]))] * -1,
                       linetype="dashed", colour="grey") +
  geom_tiplab(label=new.tiplabels.dated, fontface=tiplabel.face.datedtree, size=3, offset=1) +
  geom_tree() +
  geom_tippoint(aes(colour=lifestyle), size=2.5, show.legend=FALSE) +
  scale_colour_manual(values=col.df$colour[na.omit(match(sort(unique(metadata$lifestyle[match(species.tree$tip.label, metadata$name)])), col.df$lifestyle))],
                      na.translate=FALSE)

tiff(file=paste0("dated-tree-", Sys.Date(), ".tiff"),
     height=8, width=12, unit="in", res=300, compression="lzw")
gg.datedtree
dev.off()

tiff(file=paste0("dated-tree2-", Sys.Date(), ".tiff"),
     height=8, width=12, unit="in", res=300, compression="lzw")
gg.datedtree +
  geom_nodepoint(aes(subset=(node %in% c(64, 66))),
               size=3, colour="black") +
  geom_nodelab(aes(subset=(node %in% c(64, 66))),
               label=c("O'Donnell et\nal. 2020", "Lombard et\nal. 2015"),
               hjust=1.3,
               vjust=-0.8,
               size=3,
               colour="black")
dev.off()





gg.lifestyles <- gg.tree +
  xlim(0,0.5) +
  geom_tiplab(label=new.tiplabels, offset=0.05, fontface=tiplabel.face.speciestree, align=TRUE, linetype=NA, size=3) +
  geom_tippoint(aes(colour=lifestyle), size=3, show.legend=FALSE, shape="square") +
  scale_colour_manual(values=col.df$colour[na.omit(match(sort(unique(metadata$lifestyle[match(species.tree$tip.label, metadata$name)])), col.df$lifestyle))],
                      na.translate=FALSE)

lifestyle.grid <- metadata[c(12:19)]
rownames(lifestyle.grid) <- metadata$name

##Fix colours

col.df$lifestyle <- factor(col.df$lifestyle, levels=c("plant associated", "endophyte", "plant pathogen", "saprotroph", "insect associated", "mycoparasite", "human pathogen", "coral associated"))

gg.lifestyles.grid <- gheatmap(gg.lifestyles,
                               lifestyle.grid,
                               offset=0.005,
                               width=0.15,
                               font.size=3, 
                               colnames=FALSE,
                               hjust=0,
                               color="grey") +
  scale_fill_manual(values=col.df$colour[match(levels(col.df$lifestyle), col.df$lifestyle)],
                    breaks=levels(col.df$lifestyle),
                    labels=str_to_sentence(levels(col.df$lifestyle)),
                    na.translate=FALSE,
                    name=NULL) +
  theme(legend.position=c(0.2, 0.75),
        legend.text=element_text(size=14),
        legend.key=element_rect(size=4))

tiff(file=paste0("lifestyles-tree-", Sys.Date(), ".tiff"),
     height=8, width=12, unit="in", res=300, compression="lzw")
gg.lifestyles.grid
dev.off()


legend.tmp <- data.frame(col.df, x=1:8)

tiff(file=paste0("lifestyle-legend-", Sys.Date(), ".tiff"),
     height=3, width=5, unit="in", res=300, compression="lzw")
ggplot(legend.tmp, aes(x=x, y=1, colour=lifestyle)) +
  geom_point(size=5) +
  scale_colour_manual(values=col.df$colour[match(levels(col.df$lifestyle), col.df$lifestyle)],
                      breaks=levels(col.df$lifestyle),
                      labels=str_to_sentence(levels(col.df$lifestyle)),
                      name="Lifestyle",
                      guide=guide_legend(direction="horizontal",
                                         title.position="top",
                                         ncol=2,
                                         title.hjust=0)) +
  
  theme_minimal() +
  theme(legend.title=element_text(face="bold", size=14),
        legend.text=element_text(size=14))
dev.off()



##LIFESTYLE TEST RESULTS

#For effectors and orthogroups...
for (i in c("effectors", "orthogroups")){
  
  print(i)
  #Read in test results
  pca.result <- read.csv(paste0("lifestyle_test/", i, "/metadata.csv"))
  lifestyle.data <- read.csv(paste0("lifestyle_test/", i, "/data.csv"), row.names='genome')
  #Make distance matrix
  dist <- vegdist(lifestyle.data, method='jaccard')
  #Do permanova
  permanova <- adonis2(formula=dist ~ PC1 + PC2 + lifestyle, data=pca.result, permutations=9999)
  
  print(paste0("Phylogeny: ", round(sum(permanova$R2[1:2]) * 100), "%"))
  print(paste0("Lifestyle: ", round(sum(permanova$R2[3]) * 100), "%"))
  
  assign(paste0("permanova.", i), permanova)
}


adonis2(formula=vegdist(dist, method='jaccard'), ~ ingroup1 + ingroup2, data=orthogroups.copies.ingroup2, permutations=9999)

transposed <- as.data.frame(t(as.matrix(orthogroups.copies.ingroup2)))

anova(betadisper(vegdist(transposed, method="bray"), metadata$ingroup[match(rownames(transposed), metadata$file2)]))
anosim(transposed, metadata$ingroup[match(rownames(transposed), metadata$file2)], distance="bray", permutations=999)

#Function to perform NMDS for 1-10 dimensions to pick optimal number of dimensions
NMDS.scree <- function(x) {
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform=F, k=1)$stress), xlim=c(1, 10),ylim=c(0, 0.30), xlab="# of Dimensions", ylab="Stress")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform=F, k=i + 1)$stress))
  }
}

transposed <- transposed[,which(colSums(transposed) > 1)]
#NMDS.scree(transposed)
NMDS <- metaMDS(transposed, k=4, trymax=1000, trace=F)
#Plot stressplot
stressplot(NMDS)

metadata.nmds <- as.data.frame(scores(NMDS))
metadata.nmds$taxon <- metadata$name[match(rownames(metadata.nmds), metadata$file2)]
#Add Musa species
metadata.nmds$group <- metadata$ingroup[match(metadata.nmds$taxon, metadata$name)]

ggplot() +
  stat_ellipse(data=metadata.nmds,
               aes(x=NMDS1, y=NMDS2, fill=group, colour=group),
               alpha=0.15,
               lty="dotted",
               geom="polygon",
               type='t',
               size=0.4,
               show.legend=FALSE) +
  geom_point(data=metadata.nmds,
             aes(x=NMDS1, y=NMDS2, colour=group),
             size=2) +
  geom_label_repel(data=metadata.nmds,
                   aes(x=NMDS1, y=NMDS2, label=taxon, colour=group)) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() +
  theme(legend.text.align=0,
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        aspect.ratio=1,
        panel.grid=element_blank(),
        legend.position="top",
        legend.box="vertical",
        plot.margin=margin(0, 5, 0, 5),
        legend.margin=margin(0, 0, 0, 0),
        legend.text=element_text(size=8),
        legend.title=element_text(size=9),
        plot.title.position="plot") +
  coord_fixed()


##UPSET PLOT

for (ingroup in c(1, 2)) {
  
  if (ingroup == 1) {
    
    tree <- drop.tip(species.tree, metadata$name[metadata$ingroup != ingroup])
    
  } else {
    
    tree <- drop.tip(species.tree, outgroup)
    
  }
  
  #Plot tree
  gg.tree.upset <- ggtree(tree, branch.length="none") %<+% metadata
  
  #Make dataframe of species complex nodes
  sc.df.tmp <- data.frame(sc=names(table(metadata$speciescomplex.abb[match(tree$tip.label, metadata$name)]))[table(metadata$speciescomplex.abb[match(tree$tip.label, metadata$name)]) > 1],
                      node=NA)
  
  #Get nodes for each species complex
  for (i in 1:length(sc.df.tmp$sc)) {
    sc.df.tmp$node[i] <- getMRCA(tree, metadata$name[metadata$speciescomplex.abb == sc.df.tmp$sc[i]])
  }
  
  sc.df2.tmp <- data.frame(sc=unique(metadata$speciescomplex.abb)[is.na(match(unique(metadata$speciescomplex.abb), sc.df.tmp$sc))],
                       node=NA)
  
  gg.tree.upset.data <- gg.tree.upset[["data"]] %>%
    arrange(y)
  
  sc.df2.tmp$node <- gg.tree.upset.data$node[match(sc.df2.tmp$sc, gg.tree.upset.data$speciescomplex.abb)]
  
  sc.df.tmp <- rbind(sc.df.tmp, sc.df2.tmp)
  sc.df.tmp <- sc.df.tmp[match(na.omit(unique(gg.tree.upset.data$speciescomplex.abb)), sc.df.tmp$sc),]
  sc.df.tmp$box <- rep(c(0,1), length.out=length(sc.df.tmp$sc))
    
  for (i in 1:length(sc.df.tmp$sc)) {
  
    if (sc.df.tmp$box[i] == 0) {
      
      gg.tree.upset <- gg.tree.upset +
        geom_highlight(node=sc.df.tmp$node[sc.df.tmp$sc == sc.df.tmp$sc[i]], fill="#000000", alpha=0.05, extend=-1)
      
    }
    
    if (sc.df.tmp$node[i] > length(tree$tip.label)) {
      
      gg.tree.upset <- gg.tree.upset +
        geom_cladelabel(node=sc.df.tmp$node[sc.df.tmp$sc == sc.df.tmp$sc[i]],
                        label=sc.df.tmp$sc[i], angle=90, hjust=0.5, offset.text=-1, offset=-18, align=TRUE)
      
    }
  }
  
  #Plot for upset plot
  gg.tree.upset <- gg.tree.upset +
    #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
    #geom_tiplab() +
    #xlim(0, 0.03) +
    geom_tree() +
    geom_tippoint(aes(colour=lifestyle), size=2) +
    scale_x_reverse(limits=c(20, 0)) +
    scale_colour_manual(values=col.df$colour[na.omit(match(sort(unique(metadata$lifestyle[match(tree$tip.label, metadata$name)])), col.df$lifestyle))],
                        labels=str_to_sentence(col.df$lifestyle[na.omit(match(sort(unique(metadata$lifestyle[match(tree$tip.label, metadata$name)])), col.df$lifestyle))]),
                        na.translate=FALSE) +
    coord_cartesian(clip="off") +
    guides(colour=guide_legend(ncol=2)) +
    theme(plot.margin=margin(0, 0, 0, 0),
          legend.title=element_blank(),
          legend.box.margin=margin(0, 0, 0, 0),
          legend.margin=margin(0, 0, 0, 0),
          legend.text=element_text(size=7),
          legend.position="top")
  
  assign(paste0("gg.tree.upset.ingroup", ingroup), gg.tree.upset)
  
  #Create vector with the order of the tree tips for plotting
  tip.order <- gg.tree.upset[["data"]][order(gg.tree.upset[["data"]]$y),]
  tip.order <- rev(tip.order$label[tip.order$isTip == "TRUE"])
  #Create vector of tip labels colours
  label.cols <- rev(metadata$lifestyle[match(tip.order, metadata$name)])
  for (i in 1:length(col.df$lifestyle)) {
    label.cols[label.cols == col.df$lifestyle[i]] <- col.df$colour[i]
  }
  label.cols[is.na(label.cols)] <- "black"
  #Create vector of fontface for new genomes in this study
  label.face <- rev(metadata$new[match(tip.order, metadata$name)])
  for (i in 1:length(label.face)){
    if (label.face[i] == "Y") {
      label.face[i] <- "bold.italic"
    } else {
      label.face[i] <- "italic"
    }
  }
  #Fontface vector
  tiplabel.face <- rep("italic", length(tree$tip.label))
  tiplabel.face[which(metadata$new[match(tree$tip.label, metadata$name)] == "Y")] <- "bold.italic"
  
  #Convert effector count dataframe from binary to logical
  effector.count.logical <- as.data.frame(lapply(get(paste0("effector.count.ingroup", ingroup)), as.logical))
  
  colnames(effector.count.logical) <- metadata$name[match(colnames(effector.count.logical), metadata$file)]
  
  #Format effector count dataframe for upset plot
  upset.effectors.df <- t(effector.count.logical) %>%
    as_tibble(rownames="Taxon") %>%
    gather(Orthogroup, Member, -Taxon) %>%
    filter(Member) %>%
    select(- Member) %>%
    group_by(Orthogroup) %>%
    summarize(Taxon=list(Taxon))
  
  #Plot
  gg.upset.effectors <- ggplot(upset.effectors.df, aes(x=Taxon)) +
    geom_bar(stat="count") +
    geom_text(stat="count", aes(label=after_stat(count)), size=2, vjust=-1) +
    scale_x_upset(n_intersections=40, sets=tip.order) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.15))) +
    labs(y="Common CSEPs",
         subtitle=bquote(bold(PERMANOVA)~Phylogeny:.(round(sum(permanova.effectors$R2[1:2]) * 100))*"%"~Lifestyle:.(round(sum(permanova.effectors$R2[3]) * 100))*"%")) +
    theme_classic() +
    theme_combmatrix(combmatrix.panel.line.size=1,
                     combmatrix.label.text=element_text(face=label.face, colour=label.cols)) +
    theme(axis.title.x=element_blank(),
          plot.margin=margin(5, 0, 0, 0),
          axis.title.y=element_text(vjust=-115, size=7),
          plot.subtitle=element_text(vjust=-10, hjust=1, size=8, margin=margin(0, 0, 0, 0)))
  
  assign(paste0("gg.upset.effectors.ingroup", ingroup), gg.upset.effectors)
  
  #Convert orthogroup count dataframe from binary to logical
  orthogroups.copies.logical <- as.data.frame(lapply(get(paste0("orthogroups.copies.ingroup", ingroup)), as.logical))
  
  colnames(orthogroups.copies.logical) <- metadata$name[match(colnames(orthogroups.copies.logical), metadata$file)]
  
  #Format orthogroup count dataframe for upset plot
  upset.orthogroups.df <- t(orthogroups.copies.logical) %>%
    as_tibble(rownames="Taxon") %>%
    gather(Orthogroup, Member, -Taxon) %>%
    filter(Member) %>%
    select(- Member) %>%
    group_by(Orthogroup) %>%
    summarize(Taxon=list(Taxon))
  
  #Plot
  gg.upset.orthogroups <- ggplot(upset.orthogroups.df, aes(x=Taxon)) +
    geom_bar(stat="count") +
    geom_text(stat="count", aes(label=after_stat(count)), size=2, nudge_y=500, angle=45) +
    scale_x_upset(n_intersections=40, sets=tip.order) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.12))) +
    labs(y="Common orthogroups",
         subtitle=bquote(bold(PERMANOVA)~Phylogeny:.(round(sum(permanova.orthogroups$R2[1:2]) * 100))*"%"~Lifestyle:.(round(sum(permanova.orthogroups$R2[3]) * 100))*"%")) +
    theme_classic() +
    theme_combmatrix(combmatrix.panel.line.size=1,
                     combmatrix.label.text=element_text(face=label.face, colour=label.cols)) +
    theme(axis.title.x=element_blank(),
          plot.margin=margin(5, 0, 0, 0),
          axis.title.y=element_text(vjust=-105, size=7),
          plot.subtitle=element_text(vjust=-10, hjust=1, size=8, margin=margin(0, 0, 0, 0)))
  
  assign(paste0("gg.upset.orthogroups.ingroup", ingroup), gg.upset.orthogroups)
  
  
  ##CORE ACCESSORY SPECIFIC BARGRAPH
  
  orthogroups.copies <- get(paste0("orthogroups.copies.ingroup", ingroup))
  orthogroups.stats <- get(paste0("orthogroups.stats.ingroup", ingroup))
  
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
  secretome.df$taxon <- factor(secretome.df$taxon, levels=rev(tip.order))
  #secretome.df <- melt(secretome.df)
  
  gg.secretome.orthogroups <- ggplot(secretome.df, aes(y=taxon, x=orthogroups, fill=secretome)) +
    geom_bar(stat="identity", colour="black", size=0.5, width=0.6) +
    scale_x_reverse(expand=c(0, 0),
                    position="top") +
    scale_fill_manual(values=c("dimgrey", "darkgrey", "lightgrey"),
                      labels=c("Specific", "Accessory", "Core")) +
    theme_minimal() +
    theme(axis.title=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.y=element_blank(),
          panel.grid.major.y=element_blank(),
          legend.position="top",
          legend.text=element_text(size=10),
          legend.key.size=unit(0.3, "cm"),
          legend.spacing.x=unit(0.4, "cm"),
          legend.title=element_blank())
  
  gg.secretome.effectors <- ggplot(secretome.df, aes(y=taxon, x=effectors, fill=secretome)) +
    geom_bar(stat="identity", colour="black", size=0.5, width=0.6) +
    scale_x_reverse(expand=c(0, 0),
                    position="top") +
    scale_fill_manual(values=c("dimgrey", "darkgrey", "lightgrey"),
                      labels=c("Specific", "Accessory", "Core")) +
    theme_minimal() +
    theme(axis.title=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.y=element_blank(),
          panel.grid.major.y=element_blank(),
          legend.position="top",
          legend.text=element_text(size=10),
          legend.key.size=unit(0.3, "cm"),
          legend.spacing.x=unit(0.4, "cm"),
          legend.title=element_blank())
  
  assign(paste0("gg.secretome.orthogroups.ingroup", ingroup), gg.secretome.orthogroups)
  assign(paste0("gg.secretome.effectors.ingroup", ingroup), gg.secretome.effectors)
  assign(paste0("secretome.df.ingroup", ingroup), secretome.df)
  
}

tiff(file=paste0("upset-plot-effectors-", Sys.Date(), ".tiff"),
     width=15, height=12, units="in", res=300, compression="lzw")
cowplot::plot_grid(
  cowplot::plot_grid(NULL, gg.secretome.effectors.ingroup2 + theme(plot.margin=unit(c(1, 2, 1.5, 1), "pt")),
                     ncol=1, rel_heights=c(1, 8.2)),
  gg.upset.effectors.ingroup2,
  cowplot::plot_grid(NULL, gg.tree.upset.ingroup2 + theme(plot.margin=unit(c(1, 2, 1.5, 1), "pt")),
                     ncol=1, rel_heights=c(1, 12)),
  nrow=1, rel_widths=c(2, 3, 1))
dev.off()

tiff(file=paste0("upset-plot-orthogroups-", Sys.Date(), ".tiff"),
     width=15, height=12, units="in", res=300, compression="lzw")
cowplot::plot_grid(
  cowplot::plot_grid(NULL, gg.secretome.orthogroups.ingroup2 + theme(plot.margin=unit(c(1, 2, 1.5, 1), "pt")),
                     ncol=1, rel_heights=c(1, 8.2)),
  gg.upset.orthogroups.ingroup2,
  cowplot::plot_grid(NULL, gg.tree.upset.ingroup2 + theme(plot.margin=unit(c(1, 2, 1.5, 1), "pt")),
                     ncol=1, rel_heights=c(1, 12)),
  nrow=1, rel_widths=c(2, 3, 1))
dev.off()


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

tiff(file=paste0("pvalue-grid-", Sys.Date(), ".tiff"),
     width=6, height=3, units="in", res=300, compression="lzw")
gg.pwperm
dev.off()

tiff(file=paste0("upset-plot-effectors1-", Sys.Date(), ".tiff"),
     width=15, height=10, units="in", res=300, compression="lzw")
cowplot::plot_grid(
  cowplot::plot_grid(NULL, gg.secretome.effectors.ingroup1 + theme(plot.margin=unit(c(1, 2, 1.5, 1), "pt")),
                     ncol=1, rel_heights=c(1, 8)),
  gg.upset.effectors.ingroup1,
  cowplot::plot_grid(NULL, gg.tree.upset.ingroup1 + theme(plot.margin=unit(c(1, 2, 1.5, 1), "pt")),
                     ncol=1, rel_heights=c(1, 12)),
  nrow=1, rel_widths=c(2, 3, 1))
dev.off()

tiff(file=paste0("upset-plot-orthogroups1-", Sys.Date(), ".tiff"),
     width=15, height=10, units="in", res=300, compression="lzw")
cowplot::plot_grid(
  cowplot::plot_grid(NULL, gg.secretome.orthogroups.ingroup1 + theme(plot.margin=unit(c(1, 2, 1.5, 1), "pt")),
                     ncol=1, rel_heights=c(1, 8)),
  gg.upset.orthogroups.ingroup1,
  cowplot::plot_grid(NULL, gg.tree.upset.ingroup1 + theme(plot.margin=unit(c(1, 2, 1.5, 1), "pt")),
                     ncol=1, rel_heights=c(1, 12)),
  nrow=1, rel_widths=c(2, 3, 1))
dev.off()


tiff(file=paste0("upset-plot-", Sys.Date(), ".tiff"),
     width=10, height=12, units="in", res=300, compression="lzw")

egg::ggarrange(
  cowplot::plot_grid(
    gg.upset.orthogroups,
    cowplot::plot_grid(NULL, gg1 + theme(plot.margin=unit(c(1, -5, -5, 1), "pt")), ncol=1, rel_heights=c(1, 7.5)),
    nrow=1, rel_widths=c(3, 1)),
  cowplot::plot_grid(
    gg.upset.effectors,
    cowplot::plot_grid(NULL, gg1 + theme(plot.margin=unit(c(1, -5, -5, 1), "pt")), ncol=1, rel_heights=c(1, 7.5)),
    nrow=1, rel_widths=c(3, 1)),
  ncol=1,
  labels=c("A", "B"), label.args=list(gp=grid::gpar(font=2, cex=1.2))
)

dev.off()


##SPECIES SPECIFIC PLOT


effector.count.specific <- effector.count[!is.na(match(rownames(effector.count), orthogroups.stats$orthogroup[orthogroups.stats$secretome == "specific"])),]
effector.count.specific <- effector.count.specific[-which(colnames(effector.count.specific) == "Ilysp1_GeneCatalog_proteins_20121116")]

specific.df <- data.frame(species=colnames(effector.count.specific),
                          num=NA)

for (i in colnames(effector.count.specific)) {
  specific.df$num[specific.df$species == i] <- length(which(effector.count.specific[,i] > 0))
}

specific.df$lifestyle <- metadata$lifestyle[match(specific.df$species, metadata$file2)]
specific.df$lifestyle <- sub("-", " ", specific.df$lifestyle)

#Tukey significance testing
tukey <- TukeyHSD(aov(lm(num ~ lifestyle, data=specific.df)))
#Make dataframe for ggplot with tukey groups
tukey.df <- data.frame(multcompLetters(tukey[["lifestyle"]][,4])["Letters"])
tukey.df <- data.frame(Treatment=rownames(tukey.df), Letters=tukey.df$Letters)

ggplot(specific.df, aes(x=lifestyle, y=num, fill=lifestyle)) +
  geom_boxplot() +
  geom_text(data=tukey.df,
            aes(x=Treatment, y=Inf, label=Letters),
            family="mono",
            hjust=0.5,
            size=3,
            inherit.aes=FALSE)



##ASR

phy.ult <- phy
#phy.ult$edge.length <- 1
phy.ult <- chronos(phy.ult)
#phy.ult <- drop.tip(phy.ult, outgroup)

models <- c("lik.full", "lik.minimal", "lik.speciation", "lik.extinction", "lik.transition", "lik.specext", "lik.spectran", "lik.exttran", "fit.full", "fit.minimal", "fit.speciation", "fit.extinction", "fit.transition", "fit.specext", "fit.spectran", "fit.exttran")

models <- models[grep("lik", models)]

for (i in c("hyp1", "hyp2", "hyp3")) {
  
  print(i)
  
  states <- metadata[match(phy.ult$tip.label, metadata$name), paste0("lifestyle.", i)]
  states <- as.numeric(factor(states))
  names(states) <- metadata$name[match(phy.ult$tip.label, metadata$name)]
  
  print("Fitting ASR models")
  #p <- starting.point.musse(phy.ult, 4)
  #All varying (most complex)
  lik.full <- make.musse(phy.ult, states, 4)
  #fit.full <- find.mle(lik.full, x.init=p[argnames(lik.full)], control=list(maxit=1000000))
  #All constant (least complex)
  lik.minimal <- constrain(lik.full,
                           lambda1 ~ lambda2, lambda1 ~ lambda3, lambda1 ~ lambda4,
                           mu1 ~ mu2, mu1 ~ mu3, mu1 ~ mu4,
                           q12 ~ q21, q13 ~ q31, q14 ~ q41)
  #fit.minimal <- find.mle(lik.minimal, p[argnames(lik.minimal)], control=list(maxit=1000000))
  #Varying speciation
  lik.speciation <- constrain(lik.full,
                              mu1 ~ mu2, mu1 ~ mu3, mu1 ~ mu4,
                              q12 ~ q21, q13 ~ q31, q14 ~ q41)
  #fit.speciation <- find.mle(lik.speciation, p[argnames(lik.speciation)], control=list(maxit=1000000))
  #Varying extinction
  lik.extinction <- constrain(lik.full, 
                              lambda1 ~ lambda2, lambda1 ~ lambda3, lambda1 ~ lambda4,
                              q12 ~ q21, q13 ~ q31, q14 ~ q41)
  #fit.extinction <- find.mle(lik.extinction, p[argnames(lik.extinction)], control=list(maxit=1000000))
  #Varying transition
  lik.transition <- constrain(lik.full,
                              lambda1 ~ lambda2, lambda1 ~ lambda3, lambda1 ~ lambda4,
                              mu1 ~ mu2, mu1 ~ mu3, mu1 ~ mu4)
  #fit.transition <- find.mle(lik.transition, p[argnames(lik.transition)], control=list(maxit=1000000))
  #Varying speciation and extinction
  lik.specext <- constrain(lik.full,
                           q12 ~ q21, q13 ~ q31, q14 ~ q41)
  #fit.specext <- find.mle(lik.specext, p[argnames(lik.specext)], control=list(maxit=1000000))
  #Varying speciation and transition
  lik.spectran <- constrain(lik.full, 
                            mu1 ~ mu2, mu1 ~ mu3, mu1 ~ mu4)
  #fit.spectran <- find.mle(lik.spectran, p[argnames(lik.spectran)], control=list(maxit=1000000))
  #Varying extinction and transition
  lik.exttran <- constrain(lik.full, 
                           lambda1 ~ lambda2, lambda1 ~ lambda3, lambda1 ~ lambda4)
  #fit.exttran <- find.mle(lik.exttran, p[argnames(lik.exttran)], control=list(maxit=1000000))
  
  #results <- anova(fit.minimal, full=fit.full, speciation=fit.speciation, extinction=fit.extinction, transition=fit.transition, specext=fit.specext, spectran=fit.spectran, exttran=fit.exttran) #best model lowest AIC
  #assign(paste0("results.", i), results)
  
  for (j in models) {
    assign(paste0(j, ".", i), get(j))
  }
}

#Selected model (lowest AIC across all hypotheses)
best.model <- rownames(results)[which.min(results.hyp1$AIC + results.hyp2$AIC + results.hyp3$AIC)]
print(paste0("Selected model: ", best.model))

for (i in c("hyp1", "hyp2", "hyp3")) {
  print("Preparing pies")
  st <- asr.marginal(get(paste0("lik.", best.model, ".", i)), coef(get(paste0("fit.", best.model, ".", i))))
  pies <- as.data.frame(t(st))
  pies["node"] <- seq(1+Ntip(phy.ult), length(pies$V1)+Ntip(phy.ult), 1)
  pie.colours <- c("#009E73", "#56B4E9", "dimgrey", "#0072B2")
  plotpies <- nodepie(data=pies, color=pie.colours, cols=1:4)
  assign(paste0("plotpies.", i), plotpies)
}

gg.asr <- gg %<+% metadata +
  xlim(0, 17) +
  #FIESC
  geom_highlight(node=39, fill="#000000", alpha=0.05, extend=3.7) +
  geom_cladelabel(node=39, label="FIESC", angle=90, hjust=0.5, offset.text=0.2, offset=3.7) +
  #FSAMSC
  geom_cladelabel(node=33, label="FSAMSC", angle=90, hjust=0.5, offset.text=0.2, offset=5) +
  #FFSC
  geom_highlight(node=45, fill="#000000", alpha=0.05, extend=4.9) +
  geom_cladelabel(node=45, label="FFSC", angle=90, hjust=0.5, offset.text=0.2, offset=4.9) +
  #FOSC
  geom_cladelabel(node=29, label="FOSC", angle=90, hjust=0.5, offset.text=0.2, offset=6.3) +
  #FSSC
  geom_highlight(node=40, fill="#000000", alpha=0.05, extend=3.8) +
  geom_cladelabel(node=40, label="FSSC", angle=90, hjust=0.5, offset.text=0.2, offset=3.8) +
  geom_tree(size=2) +
  geom_tiplab(fontface=tiplabel.face, size=5, offset=0.2, colour="black") +
  scale_colour_manual(values=c("#009E73", "#56B4E9", "dimgrey", "#0072B2"), labels=c("Endophyte", "Insect-mutualist", "Plant pathogen", "Saprotroph"), na.translate=FALSE, guide=guide_legend(title=NULL)) +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=15))

gg.asr.hyp1 <- inset(gg.asr, plotpies.hyp1, width=0.09, height=0.09)
gg.asr.hyp1 <- gg.asr.hyp1 +
  geom_tippoint(aes(colour=lifestyle.hyp1), size=5)
gg.asr.hyp3 <- inset(gg.asr, plotpies.hyp3, width=0.09, height=0.09)
gg.asr.hyp3 <- gg.asr.hyp3 +
  geom_tippoint(aes(colour=lifestyle.hyp3), size=5)

gg.asr.all <- ggarrange(gg.asr.hyp1, gg.asr.hyp3,
                        nrow=1,
                        labels=c("Hypothesis 1", "Hypothesis 3"),
                        hjust=0,
                        font.label=list(size=8, face="bold"),
                        common.legend=TRUE)
#annotate_figure(gg.all, fig.lab="True sampling proportions estimated", fig.lab.face="bold")

#tiff(file=paste0("ASRtree-", Sys.Date(), ".tiff"), height=7, width=12, unit="in", res=300)
gg.asr.all
#dev.off()




##SELECTION

core.SC.orthogroups <- Reduce(intersect,
                              list(orthogroups.stats.ingroup0$orthogroup[which(orthogroups.stats.ingroup0$copy_number == "single")],
                                   orthogroups.stats.ingroup0$orthogroup[which(orthogroups.stats.ingroup0$secretome == "core")]))
  
absrel.p <- list()

for (i in core.SC.orthogroups) {
  
  absrel.results <- tryCatch(fromJSON(paste0("selection/hyphy/absrel/", i, "_aBSREL.json")), error=function(e) NULL)
  
  if (!is.null(absrel.results)) {
    #assign(paste0("absrel.results.", i), absrel.results)
    
    for (j in 1:length(names(absrel.results[["branch attributes"]][["0"]]))) {
      absrel.p[[names(absrel.results[["branch attributes"]][["0"]][j])]][i] <- absrel.results[["branch attributes"]][["0"]][[names(absrel.results[["branch attributes"]][["0"]])[j]]][["Corrected P-value"]]
    }
  }
}


blah <- read.tree(text=paste0(absrel.results.OG0005792$input$trees[[1]], ";"))
cowplot::plot_grid(
  ggtree(blah) +
    xlim(c(0,20)) + 
    geom_tiplab() +
    geom_nodelab(),
  ggtree(species.tree, branch.length="none") +
    xlim(c(0,25)) + 
    geom_tiplab() +
    geom_text2(aes(subset=!isTip, label=node), hjust=-0.1)
)

core.SC.mixed <- Reduce(intersect,
                        list(orthogroups.stats.ingroup0$orthogroup[which(orthogroups.stats.ingroup0$copy_number == "single")],
                             orthogroups.stats.ingroup0$orthogroup[which(orthogroups.stats.ingroup0$secretome == "core")],
                             orthogroups.stats.ingroup0$orthogroup[which(orthogroups.stats.ingroup0$effector != "")]))

absrel.tree <- root(iqtree, "Ilysp1_GeneCatalog_proteins_20121116", resolve.root=TRUE, edgelabel=TRUE)

absrel.df <- read.csv("selection/hyphy/absrel/absrel_nodes.csv")
absrel.df <- rbind(absrel.df, data.frame(node=1:length(absrel.tree$tip.label), absrel=absrel.tree$tip.label))

absrel.df$absrel[which(absrel.df$absrel %in% metadata$name)] <- metadata$tip[match(absrel.df$absrel[which(absrel.df$absrel %in% metadata$name)], metadata$name)]

absrel.p.count <- lengths(lapply(absrel.p, function(x) which(x <=0.05)))
absrel.df$num <- absrel.p.count[match(absrel.df$absrel, names(absrel.p.count))]

absrel.p.names <- unlist(lapply(absrel.p, function(x) paste(names(which(x <=0.05)), collapse=" ")))
absrel.p.names <- absrel.p.names[absrel.p.names != ""]
#absrel.p.names <- gsub("OG000", "", absrel.p.names)

absrel.df$orthos <- NA
absrel.df$orthos[match(names(absrel.p.names[grep(core.SC.mixed[which(core.SC.mixed %in% absrel.p.names)], absrel.p.names)]), absrel.df$absrel)] <- gsub("OG000", "", core.SC.mixed[which(core.SC.mixed %in% absrel.p.names)])

#absrel.df$orthos <- absrel.p.names[match(absrel.df$absrel, names(absrel.p.names))]
#absrel.df$orthos[absrel.df$orthos == ""] <- NA
#absrel.df$orthos <- gsub('(?=(?:.{10})+$)', "\n", absrel.df$orthos, perl = TRUE)

species.tree <- root(iqtree, "Ilysp1_GeneCatalog_proteins_20121116", resolve.root=TRUE, edgelabel=TRUE)
species.tree$tip.label <- metadata$name[match(species.tree$tip.label, metadata$tip)]

length(which(as.numeric(species.tree$node.label) > 95)) / length(na.omit(as.numeric(species.tree$node.label)))

gg.selection <- ggtree(species.tree, branch.length="none") %<+% metadata

#Create vector with the order of the tree tips for plotting
gg.tree.data.iq <- gg.selection[["data"]] %>%
  arrange(y)
tip.order.speciestree <- rev(gg.tree.data.iq$label[gg.tree.data.iq$isTip == "TRUE"])
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
tiplabel.face.speciestree <- rep("italic", length(species.tree$tip.label))
tiplabel.face.speciestree[which(metadata$new[match(species.tree$tip.label, metadata$name)] == "Y")] <- "bold.italic"

new.tiplabels.iq <- species.tree$tip.label

for (i in 1:length(new.tiplabels.iq)) {
  if (metadata$old.name[match(new.tiplabels.iq, metadata$name)][i] != "") {
    new.tiplabels.iq[i] <- paste0(new.tiplabels.iq[i],
                                  " (=", metadata$old.name[match(new.tiplabels.iq, metadata$name)][i], ")")
  }
}

new.tiplabels.iq <- sub("Fusarium", "F.", new.tiplabels.iq)

#Make dataframe of species complex nodes
sc.df.iq <- data.frame(sc=names(table(metadata$speciescomplex.abb[match(species.tree$tip.label, metadata$name)]))[table(metadata$speciescomplex.abb[match(species.tree$tip.label, metadata$name)]) > 1],
                       node=NA)

#Get nodes for each species complex
for (i in 1:length(sc.df.iq$sc)) {
  sc.df.iq$node[i] <- getMRCA(species.tree, metadata$name[metadata$speciescomplex.abb == sc.df.iq$sc[i]])
}

sc.df.iq2 <- data.frame(sc=unique(metadata$speciescomplex.abb)[is.na(match(unique(metadata$speciescomplex.abb), sc.df.iq$sc))],
                     node=NA)

sc.df.iq2$node <- gg.tree.data.iq$node[match(sc.df.iq2$sc, gg.tree.data.iq$speciescomplex.abb)]

sc.df.iq <- rbind(sc.df.iq, sc.df.iq2)
sc.df.iq <- sc.df.iq[match(na.omit(unique(gg.tree.data.iq$speciescomplex.abb)), sc.df.iq$sc),]
sc.df.iq$box <- rep(c(0,1), length.out=length(sc.df.iq$sc))

for (i in 1:length(sc.df.iq$sc)) {
  
  if (sc.df.iq$box[i] == 1) {
    
    gg.selection <- gg.selection +
      geom_highlight(node=sc.df.iq$node[sc.df.iq$sc == sc.df.iq$sc[i]], fill="#000000", alpha=0.05, extend=10)
    
  } else {
    
    gg.selection <- gg.selection +
      geom_highlight(node=sc.df.iq$node[sc.df.iq$sc == sc.df.iq$sc[i]], fill="#000000", alpha=0.1, extend=10)
    
  }
  
}

gg.selection <- gg.selection %<+% absrel.df +
  geom_tree(size=2, aes(colour=num)) +
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
  geom_tippoint(aes(colour=lifestyle), size=4) +
  scale_colour_manual(values=col.df$colour[na.omit(match(sort(unique(metadata$lifestyle[match(species.tree$tip.label, metadata$name)])), col.df$lifestyle))],
                      labels=str_to_sentence(col.df$lifestyle[na.omit(match(sort(unique(metadata$lifestyle[match(species.tree$tip.label, metadata$name)])), col.df$lifestyle))]),
                      na.translate=FALSE,
                      guide=guide_legend(title="Lifestyle",
                                         ncol=2,
                                         title.hjust=0)) +
  new_scale_colour() +
  #geom_tiplab(label=new.tiplabels, fontface=tiplabel.face.speciestree, size=4, offset=0.1, colour="black") +
  geom_label2(aes(x=branch, label=num, fill=num), size=2, label.size=0, fontface="bold") +
  #geom_label2(aes(x=branch, label=orthos, fill=num), size=2, label.size=0, fontface="bold") +
  scale_fill_gradient2(trans="pseudo_log",
                         low="black", mid="#F0E442", high="#CC79A7", midpoint=1,
                         breaks=c(0, 10, 100),
                         limits=c(0, 100),
                         labels=c(0 , 10, 100),
                         na.value="dimgrey",
                         guide=guide_colourbar(title="Orthogroups showing\npositive selection",
                                               title.position="top",
                                               direction="horizontal")) +
  theme(legend.box="vertical",
        #legend.position=c(0.15, 0.85),
        legend.title=element_text(size=14, face="bold"),
        legend.text=element_text(size=14))

tiff(file=paste0("presentationtree1-", Sys.Date(), ".tiff"),
     height=10, width=14, unit="in", res=300, compression="lzw")
gg.selection +
  geom_tiplab(label=new.tiplabels.iq, fontface=tiplabel.face.speciestree, size=4, offset=0.1, colour="black") +
  new_scale_fill() +
  coord_cartesian(clip="off") +
  theme(legend.position=c(0.15, 0.85),
        legend.text=element_text(size=14),
        legend.key=element_rect(size=4),
        plot.margin=margin(0, 5.5, 0, 50, unit="pt"))
dev.off()

tiff(file=paste0("presentationtree2-", Sys.Date(), ".tiff"),
     height=10, width=14, unit="in", res=300, compression="lzw")
gg.selection +
  geom_tiplab(label=new.tiplabels.iq, fontface=tiplabel.face.speciestree, size=4, offset=0.1, colour="black") +
  geom_nodepoint(aes(subset=(node %in% c(84, 123))),
                 size=4) +
  geom_nodelab(aes(subset=(node %in% c(84, 123))),
               label=c("Lombard et al. 2015\ngeneric concept", "O'Donnell et al. 2020\ngeneric concept"),
               hjust=1.1,
               vjust=-0.5,
               size=4) +
  new_scale_fill() +
  coord_cartesian(clip="off") +
  theme(legend.position=c(0.15, 0.85),
        legend.text=element_text(size=14),
        legend.key=element_rect(size=4),
        plot.margin=margin(0, 5.5, 0, 50, unit="pt"))
dev.off()

tiff(file=paste0("presentationtree3-", Sys.Date(), ".tiff"),
     height=10, width=14, unit="in", res=300, compression="lzw")
gg.pres3 <- ggtree(species.tree, branch.length="none") %<+% metadata +
  xlim(0, 25) +
  geom_tippoint(aes(colour=lifestyle), size=4) +
  scale_colour_manual(values=col.df$colour[na.omit(match(sort(unique(metadata$lifestyle[match(species.tree$tip.label, metadata$name)])), col.df$lifestyle))],
                      labels=str_to_sentence(col.df$lifestyle[na.omit(match(sort(unique(metadata$lifestyle[match(species.tree$tip.label, metadata$name)])), col.df$lifestyle))]),
                      na.translate=FALSE,
                      guide=guide_legend(title="Lifestyle",
                                         ncol=2,
                                         title.hjust=0)) +
  geom_tiplab(label=new.tiplabels, fontface=tiplabel.face.speciestree, size=4, offset=2, colour="black") +
  theme(legend.position="none")
gg.pres3
dev.off()

tiff(file=paste0("presentationtree4-", Sys.Date(), ".tiff"),
     height=10, width=14, unit="in", res=300, compression="lzw")
gg.pres4 <- gg.pres3 +
  geom_tiplab(label=new.tiplabels, fontface=tiplabel.face.speciestree, size=4, offset=2, colour="black") +
  new_scale_fill() +
  coord_cartesian(clip="off")

gheatmap(gg.pres3,
         lifestyle.grid,
         offset=0.1,
         width=0.1,
         font.size=3, 
         colnames=FALSE,
         hjust=0,
         color="grey",) +
  scale_fill_manual(values=col.df$colour[match(levels(col.df$lifestyle), col.df$lifestyle)],
                    breaks=levels(col.df$lifestyle),
                    labels=str_to_sentence(levels(col.df$lifestyle)),
                    na.translate=FALSE,
                    guide=FALSE) +
  theme(legend.position="none")


dev.off()


tiff(file=paste0("postertree-", Sys.Date(), ".tiff"),
     height=10, width=17, unit="in", res=300, compression="lzw")
gg.poster <- gg.selection +
  geom_tiplab(label=new.tiplabels, fontface=tiplabel.face.speciestree, size=4, offset=2, colour="black") +
  geom_nodepoint(aes(subset=(node %in% c(84, 123))),
                            size=4) +
  geom_nodelab(aes(subset=(node %in% c(84, 123))),
               label=c("Lombard et al. 2015\ngeneric concept", "O'Donnell et al. 2020\ngeneric concept"),
               hjust=1.1,
               vjust=-0.5,
               size=4) +
  new_scale_fill() +
  coord_cartesian(clip="off")

gheatmap(gg.poster,
         lifestyle.grid,
         offset=0.1,
         width=0.1,
         font.size=3, 
         colnames=FALSE,
         hjust=0,
         color="grey",) +
  scale_fill_manual(values=col.df$colour[match(levels(col.df$lifestyle), col.df$lifestyle)],
                    breaks=levels(col.df$lifestyle),
                    labels=str_to_sentence(levels(col.df$lifestyle)),
                    na.translate=FALSE,
                    guide=FALSE) +
  theme(legend.position="none",
        legend.text=element_text(size=14),
        legend.key=element_rect(size=4),
        plot.margin=margin(0, 15, 0, 40, unit="pt"))
dev.off()

gg %<+% absrel.df +
  scale_colour_manual(values=c("#009E73", "#56B4E9", "dimgrey", "#0072B2"), 
                      labels=c("Endophyte", "Insect-mutualist", "Plant pathogen", "Saprotroph"), na.translate=FALSE,
                      guide=guide_legend(title="Lifestyle",                                                                                                       ncol=2,
                                         title.position="top",
                                         title.hjust=0)) +
  new_scale_colour() +
  aes(colour=num) +
  #geom_label2(aes(subset=(num>0), x=branch, label=num), size=3, fontface="bold") +
  scale_colour_gradient(low="#F0E442", high="#CC79A7") +
  new_scale_colour() +
  geom_tiplab(size=5, colour="black") +
  geom_label2(aes(x=branch, label=orthos, fill=num), size=2, label.size=0, fontface="bold") +
  scale_fill_gradient(low="#F0E442", high="#CC79A7",
                      guide=guide_colourbar(title="Effectors showing positive selection",
                                            title.position="top")) +
  theme(legend.box="horizontal",
        legend.title=element_text(size=10, face="bold"),
        legend.text=element_text(size=10))


##BUSTED

busted.df <- data.frame(orthogroup=core.SC.orthogroups, p=NA, effector=NA)
busted.df$effector[match(core.SC.mixed, busted.df$orthogroup)] <- "Y"

for (i in busted.df$orthogroup) {
  
  busted.results <- tryCatch(fromJSON(paste0("selection/hyphy/busted/", i, "_BUSTED.json")), error=function(e) NULL)
  
  if (!is.null(busted.results)) {
    #assign(paste0("absrel.results.", i), absrel.results)
    
    busted.df$p[match(i, busted.df$orthogroup)] <- busted.results$`test results`$`p-value`
  }
}

busted.sum <- data.frame(group=c("CSEPs", "Other"),
                         total=c(length(which(busted.df$effector == "Y")),
                                 length(which(is.na(busted.df$effector)))),
                         signif=c(length(which(busted.df$p[which(busted.df$effector == "Y")] < 0.05)),
                                  length(which(busted.df$p[which(is.na(busted.df$effector))] < 0.05))))
busted.sum$label <- paste0(round(busted.sum$signif/busted.sum$total*100), "%")
busted.sum <- melt(busted.sum)
busted.sum$label[busted.sum$variable == "total"] <- NA

busted.label <- data.frame(x=1.5, y=1.5, label="Total single-copy\ncore orthogroups",
                           group=factor("Other", levels=c("CSEPs", "Other")))

ggplot(busted.sum, aes(x=1, y=value)) +
  facet_grid(~group, switch="x") +
  geom_bar(aes(fill=variable), stat="identity") +
  scale_fill_manual(values=c("grey", "dimgrey"))
  new_scale(new_aes="size") +
  geom_text(aes(label=label, size=value), colour="white") +
  geom_text(data=busted.label, aes(x=x, y=y, label=label), inherit.aes=FALSE) +
  scale_size(range=c(2,15), guide="none") +
  scale_x_continuous(limits=c(0,3)) +
  scale_y_continuous(limits=c(0,3)) +
  #theme_void() +
  theme(legend.position="none",
        strip.text=element_text(face="bold", size=15))


##MEME

species.tree <- root(species.tree, "Ilysp1_GeneCatalog_proteins_20121116", edgelabel=TRUE)

msa.labels <- sub("[[:space:]]+\\(.*", "", new.tiplabels)

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



ggarrange(
  ggtree(meme.tree.gene) +
    xlim(c(0,17)) + 
    geom_tiplab() +
    geom_nodelab(),
  ggtree(gene.tree.5992, branch.length="none") +
    xlim(c(0,20)) + 
    geom_tiplab() +
    #geom_nodelab() +
    geom_text2(aes(subset=!isTip, label=node), hjust=-1))

ggtree(gene.tree.5992, branch.length="none") +
  xlim(c(0,20)) + 
  geom_tiplab() +
  #geom_nodelab() +
  geom_text2(aes(subset=!isTip, label=node), hjust=-1)


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
  
  assign(paste0("gg.meme.", i), gg.meme)
  
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

