setwd("D:/Documents/GitHub/FusariumEndophytes/")

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
library(diversitree)
library(colorspace)
library(jsonlite)
library(scales)
library(patchwork)
library(seqmagick)
library(ggmsa)
library(multcompView)
library(ggstance)

#Colour palette
scales::show_col(colorblind_pal()(8))

#ADD read in effector count dataframe from parser script

specific.df <- data.frame(species=rep(colnames(effector.count)[-c(1,2)], 2),
                          copy=rep(c("single", "multi"), each=length(colnames(effector.count)[-c(1,2)])),
                          num=NA)

for (i in 1:length(unique(specific.df$species))) {
  specific.df[specific.df$species == specific.df$species[i] & specific.df$copy == "single", "num"] <- length(which(effector.count.SP[,specific.df$species[i]] == 1))
  specific.df[specific.df$species == specific.df$species[i] & specific.df$copy == "multi", "num"] <- length(which(effector.count.SP[,specific.df$species[i]] > 1))
}

#Read in sample metadata
metadata <- read.csv("metadata.csv")


specific.df$lifestyle <- metadata$lifestyle.hyp1[match(specific.df$species, metadata$file)]
specific.df <- specific.df[-which(specific.df$lifestyle == ""),]

#Tukey significance testing
tukey <- TukeyHSD(aov(lm(num ~ lifestyle, data=specific.df[specific.df$copy == "single",])))
#Make dataframe for ggplot with tukey groups
tukey.df <- data.frame(multcompLetters(tukey[["num"]][,4])["Letters"])
tukey.df <- data.frame(Treatment=rownames(tukey.df), Letters=tukey.df$Letters)

ggplot(specific.df, aes(x=lifestyle, y=num, fill=lifestyle)) +
  facet_wrap(copy ~., scales="free") +
  geom_boxplot() +
  geom_text(data=tukey.df,
            aes(x=Treatment, y=Inf, label=Letters),
            family="mono",
            hjust=1,
            angle=90,
            size=3)


effector.count.SC.SP$secretome



#TANGLEGRAM

library(ape)
library(dendextend)

astral <- read.tree("phylogenomics/species_tree/astral/fus_astral_proteins_62T.tre")
astral$edge.length <- rep(1, length(astral$edge.length))
raxml <- read.tree("RAxML/RAxML_bipartitions.fus_proteins_27T_bsadd")

metadata <- read.csv("metadata.csv")

raxml$tip.label <- metadata$name[match(raxml$tip.label, metadata$tip)]
astral$tip.label <- metadata$name[match(astral$tip.label, metadata$tip)]
outgroup <- "Ilyonectria sp."

for (i in c("astral", "raxml")) {
  tree <- get(i)
  tree <- root(tree, outgroup, resolve.root=TRUE, edgelabel=TRUE)
  tree <- chronos(tree)
  tree <- as.dendrogram(tree)
  assign(paste0(i, ".dend"), tree)
}

#Untangle dendrograms
dend <- untangle(astral.dend, raxml.dend, method="step2side")

line.col <- rep("grey", length(labels(raxml.dend)))
line.col[c(59:62)] <- "black"

tiff(file=paste0("tanglegram-", Sys.Date(), ".tiff"), height=8, width=8, units="in", res=300)

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
           color_lines=line.col,
           rank_branches=TRUE)

dev.off()




##LIFESTYLE VERSUS PHYLOGENY TEST

#Read in tree
phy <- read.tree("../../Phylogenomic analysis/Proteins/RAxML/RAxML_bipartitions.fus_proteins_27T_bsadd")
#Read in tip label data
#tip.labels <- read.csv("../../Phylogenomic analysis/Proteins/tip_labels.csv")
#Replace tip labels
phy$tip.label <- metadata$name[match(phy$tip.label, metadata$tip)]
#Root tree
outgroup <- "Ilyonectria sp."
phy <- root(phy, outgroup, resolve.root=TRUE, edgelabel=TRUE)
#Remove outgroup from tree and write to file
#phy.ingroup <- drop.tip(phy, outgroup)
#write.tree(phy.ingroup, "RAxML_bipartitions.fus_proteins_27T_bsadd_rooted_ingroup")

#Replace sample IDs with names
colnames(effector.count.SC.SP) <- c("secretome", "mixed", metadata$name[match(colnames(effector.count.SC.SP)[c(3:29)], metadata$file)])
colnames(orthogroups.copies) <- metadata$name[match(colnames(orthogroups.copies), metadata$file)]

#Effectors
#Transpose count dataframe (excluding outgroup)
lifestyle.test.effectors <- as.data.frame(t(effector.count.SC.SP[c(3:23, 25:29)]))
#Add column with names
lifestyle.test.effectors$genome <- rownames(lifestyle.test.effectors)
#Add column with lifestyle
lifestyle.test.effectors$lifestyle <- gsub(" ", "", metadata$lifestyle.hyp1[match(rownames(lifestyle.test.effectors), metadata$name)])
#Replace spaces and hyphens to match tree labels
lifestyle.test.effectors$lifestyle <- gsub("-", "", lifestyle.test.effectors$lifestyle)
lifestyle.test.effectors$genome <- gsub(" ", "_", lifestyle.test.effectors$genome)
#Reorder columns
lifestyle.test.effectors <- lifestyle.test.effectors[c(402,403,1:401)]
#Write to file
#write.csv(lifestyle.test.effectors, "lifestyle-test-effectors.csv", row.names=FALSE)

#Repeat for orthogroups
lifestyle.test.orthogroups <- as.data.frame(t(orthogroups.copies[c(1:21, 23:27)]))
lifestyle.test.orthogroups$genome <- rownames(lifestyle.test.orthogroups)
lifestyle.test.orthogroups$lifestyle <- gsub(" ", "", metadata$lifestyle.hyp1[match(rownames(lifestyle.test.orthogroups), metadata$name)])
lifestyle.test.orthogroups$lifestyle <- gsub("-", "", lifestyle.test.orthogroups$lifestyle)
lifestyle.test.orthogroups$genome <- gsub(" ", "_", lifestyle.test.orthogroups$genome)
lifestyle.test.orthogroups <- lifestyle.test.orthogroups[c(23474,23475,1:23473)]
#write.csv(lifestyle.test.orthogroups, "lifestyle-test-orthogroups.csv", row.names=FALSE)

#Run python script

#For effectors and orthogroups...
for (i in c("effectors", "orthogroups")){
  
  print(i)
  #Read in test results
  pca.result <- read.csv(paste0("D:/Documents/Bio Programmes/Effect-Of-Biological-Categories-On-Genomes-Composition-master/Effect-Of-Biological-Categories-On-Genomes-Composition-master/raxml/fus_", i, "/metadata.csv"))
  lifestyle.data <- read.csv(paste0("D:/Documents/Bio Programmes/Effect-Of-Biological-Categories-On-Genomes-Composition-master/Effect-Of-Biological-Categories-On-Genomes-Composition-master/raxml/fus_", i, "/data.csv"), row.names='genome')
  #Make distance matrix
  dist <- vegdist(lifestyle.data, method='jaccard')
  #Do permanova
  permanova <- adonis2(formula=dist ~ PC1 + PC2 + lifestyle, data=pca.result, permutations=9999)
  
  print(paste0("Phylogeny: ", round(sum(permanova$R2[1:2]) * 100), "%"))
  print(paste0("Lifestyle: ", round(sum(permanova$R2[3]) * 100), "%"))
  
  assign(paste0("permanova.", i), permanova)
}

##CORE ACCESSORY SPECIFIC BARGRAPH

#Assess which orthogroups are in all species (core), one species (specific), or some species (accessory)
secretome <- vector(mode="character", length=length(rownames(orthogroups.copies)))

print("Assigning orthogroups as core, accessory or specific")
#Create bar to show progress
progress.bar <- txtProgressBar(1, length(rownames(orthogroups.copies)), initial=0, char="=", style=3)
for (j in 1:length(rownames(orthogroups.copies))) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, j)
  
  if (length(which(orthogroups.copies[j,] == 0)) == (length(colnames(orthogroups.copies)) - 1)) {
    secretome[j] <- "specific"
  }
  if (length(which(orthogroups.copies[j,] == 0)) == 0) {
    secretome[j] <- "core"
  }
  if (length(which(orthogroups.copies[j,] == 0)) < (length(colnames(orthogroups.copies)) - 1) && length(which(orthogroups.copies[j,] == 0)) > 0) {
    secretome[j] <- "accessory"
  }
}
close(progress.bar)

#Add column to effector count dataframe with whether orthogroup is SP-only or SP-mixed
orthogroups.copies$mixed <- NA

for (j in 1:length(rownames(orthogroups.copies))) {
  if(grepl("SP-mixed", mixed[j]) == TRUE) {
    orthogroups.copies$mixed[j] <- "SP-mixed"
  } else {
    orthogroups.copies$mixed[j] <- "SP-only"
  }
}

#Add column to effector count dataframe with whether orthogroup core, specific or accessory
orthogroups.copies$secretome <- secretome
#Reorder dataframe columns
orthogroups.copies <- orthogroups.copies %>% select(secretome, mixed, everything())

orthogroups.df <- data.frame(taxon=rep(colnames(orthogroups.copies)[-c(1,2)], each=3), 
                             secretome=rep(c("specific", "accessory", "core"), length(colnames(effector.count)[-c(1,2)])),
                             count=NA)

for (i in unique(orthogroups.df$taxon)) {
  for (j in unique(orthogroups.df$secretome)) {
    orthogroups.df$count[intersect(which(orthogroups.df$taxon == i), which(orthogroups.df$secretome == j))] <- table(orthogroups.copies$secretome[orthogroups.copies[i] > 0])[j]
  }
}

orthogroups.df$taxon <- metadata$name[match(orthogroups.df$taxon, metadata$file)]
orthogroups.df$secretome <- factor(orthogroups.df$secretome, levels=c("specific", "accessory", "core"))


secretome.df <- data.frame(taxon=rep(colnames(effector.count)[-c(1,2)], each=3), 
                 secretome=rep(c("specific", "accessory", "core"), length(colnames(effector.count)[-c(1,2)])),
                 count=NA)


for (i in unique(secretome.df$taxon)) {
  for (j in unique(secretome.df$secretome)) {
    secretome.df$count[intersect(which(secretome.df$taxon == i), which(secretome.df$secretome == j))] <- table(effector.count$secretome[effector.count[i] > 0])[j]
  }
}

secretome.df$taxon <- metadata$name[match(secretome.df$taxon, metadata$file)]
secretome.df$secretome <- factor(secretome.df$secretome, levels=c("specific", "accessory", "core"))

ggplot(secretome.df, aes(y=taxon, x=count, fill=secretome)) +
  geom_bar(stat="identity")

test <- ggtree(astral) + geom_tiplab()
test1 <- facet_plot(test + xlim_tree(50),
           panel="Effectors",
           data=secretome.df,
           geom=geom_barh,
           aes(x=count, fill=secretome),
           stat='identity')

test2 <- facet_plot(test1,
             panel="Orthogroups",
             data=orthogroups.df,
             geom=geom_barh,
             aes(x=count, fill=secretome),
             stat='identity') +
  theme(axis.text.x=element_text(),
        axis.line.x=element_line(),
        axis.ticks.x=element_line())


##UPSET PLOT

#Make dataframe of lifestyle colours
col.df <- data.frame(lifestyle=c("endophyte", "insect-mutualist", "plant pathogen", "saprotroph"), colour=c("#009E73", "#56B4E9", "dimgrey", "#0072B2"))

#Make dataframe of species complex nodes
sc.df <- data.frame(sc=names(table(metadata$speciescomplex.abb))[table(metadata$speciescomplex.abb) > 1], node=NA)

#Get nodes for each species complex
for (i in 1:length(sc.df$sc)) {
  sc.df$node[i] <- getMRCA(phy, metadata$name[metadata$speciescomplex.abb == sc.df$sc[i]])
}

#Plot tree
gg <- ggtree(phy, branch.length="none")

#Reorder to match BUSCO phylogeny
gg <- flip(gg, 33, 39)
gg <- flip(gg, 32, 30)

#Plot for upset plot
gg1 <- gg %<+% metadata +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  #geom_tiplab() +
  #xlim(0, 11) +
  #FIESC
  geom_highlight(node=sc.df$node[sc.df$sc == "FIESC"], fill="#000000", alpha=0.05) +
  geom_cladelabel(node=sc.df$node[sc.df$sc == "FIESC"],
                  label="FIESC", angle=90, hjust=0.5, offset.text=-0.5, offset=-1.3) +
  #FSAMSC
  geom_cladelabel(node=sc.df$node[sc.df$sc == "FSAMSC"],
                  label="FSAMSC", angle=90, hjust=0.5, offset.text=-0.5, offset=-4.36) +
  #FFSC
  geom_highlight(node=sc.df$node[sc.df$sc == "FFSC"], fill="#000000", alpha=0.05) +
  geom_cladelabel(node=sc.df$node[sc.df$sc == "FFSC"],
                  label="FFSC", angle=90, hjust=0.5, offset.text=-0.5, offset=-6.4) +
  #FOSC
  geom_cladelabel(node=sc.df$node[sc.df$sc == "FOSC"],
                  label="FOSC", angle=90, hjust=0.5, offset.text=-0.5, offset=-1.3) +
  #FSSC
  geom_highlight(node=sc.df$node[sc.df$sc == "FSSC"], fill="#000000", alpha=0.05) +
  geom_cladelabel(node=sc.df$node[sc.df$sc == "FSSC"],
                  label="FSSC", angle=90, hjust=0.5, offset.text=-0.5, offset=-5.3) +
  geom_tree() +
  geom_tippoint(aes(colour=lifestyle.hyp1), size=2) +
  scale_x_reverse(limits=c(11.5, 0)) +
  scale_colour_manual(values=c("#009E73", "#56B4E9", "dimgrey", "#0072B2"),
                      labels=c("Endophyte", "Insect-mutualist", "Plant pathogen", "Saprotroph"),
                      na.translate=FALSE) +
  guides(colour=guide_legend(ncol=2)) +
  theme(plot.margin=margin(0, 0, 5, 0),
        legend.title=element_blank(),
        legend.box.margin=margin(0, 0, 0, 0),
        legend.margin=margin(0, 0, 0, 0),
        legend.text=element_text(size=8),
        legend.position="top")

#Create vector with the order of the tree tips for plotting
tip.order <- gg[["data"]][order(gg[["data"]]$y),]
tip.order <- rev(tips$label[tips$isTip == "TRUE"])
#Create vector of tip labels colours
label.cols <- rev(metadata$lifestyle.hyp1[match(tip.order, metadata$name)])
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
tiplabel.face <- rep("italic", length(phy$tip.label))
tiplabel.face[which(metadata$new[match(phy$tip.label, metadata$name)] == "Y")] <- "bold.italic"

#Convert effector count dataframe from binary to logical
effector.count.SC.SP[3:29] <- lapply(effector.count.SC.SP[3:29], as.logical)

#Format effector count dataframe for upset plot
upset.effectors.df <- t(effector.count.SC.SP[3:29]) %>%
  as_tibble(rownames="Species") %>%
  gather(Orthogroup, Member, -Species) %>%
  filter(Member) %>%
  select(- Member) %>%
  group_by(Orthogroup) %>%
  summarize(Species=list(Species))

#Plot
gg.upset.effectors <- ggplot(upset.effectors.df, aes(x=Species)) +
  geom_bar(stat="count") +
  geom_text(stat="count", aes(label=after_stat(count)), size=2, vjust=-1) +
  scale_x_upset(n_intersections=40, sets=tip.order) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.15))) +
  labs(y="Common effector orthogroups",
       subtitle=bquote(bold(PERMANOVA)~Phylogeny:.(round(sum(permanova.effectors$R2[1:2]) * 100))*"%"~Lifestyle:.(round(sum(permanova.effectors$R2[3]) * 100))*"%")) +
  theme_classic() +
  theme_combmatrix(combmatrix.panel.line.size=1,
                   combmatrix.label.text=element_text(face=label.face, colour=label.cols)) +
  theme(axis.title.x=element_blank(),
        plot.margin=margin(5, 0, 0, 0),
        axis.title.y=element_text(vjust=-100, size=7),
        plot.subtitle=element_text(vjust=-10, hjust=1, size=8, margin=margin(0, 0, 0, 0)))

#Convert orthogroup count dataframe from binary to logical
orthogroups.copies[1:27] <- lapply(orthogroups.copies[1:27], as.logical)

#Format orthogroup count dataframe for upset plot
upset.orthogroups.df <- t(orthogroups.copies[1:27]) %>%
  as_tibble(rownames="Species") %>%
  gather(Orthogroup, Member, -Species) %>%
  filter(Member) %>%
  select(- Member) %>%
  group_by(Orthogroup) %>%
  summarize(Species=list(Species))

#Plot
gg.upset.orthogroups <- ggplot(upset.orthogroups.df, aes(x=Species)) +
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
        axis.title.y=element_text(vjust=-90, size=7),
        plot.subtitle=element_text(vjust=-10, hjust=1, size=8, margin=margin(0, 0, 0, 0)))

tiff(file=paste0("upset-plot-", Sys.Date(), ".tiff"), width=10, height=12, units="in", res=300)

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

absrel.p <- list()

for (i in rownames(effector.count.SC.SP[effector.count.SC.SP$secretome == "core",])) {
  
  absrel.results <- fromJSON(paste0("aBSREL_results/raxml/", i, "_aBSREL.json"))
  assign(paste0("absrel.results.", i), absrel.results)
  
  for (j in 1:length(names(absrel.results[["branch attributes"]][["0"]]))) {
    absrel.p[[names(absrel.results[["branch attributes"]][["0"]][j])]][i] <- absrel.results[["branch attributes"]][["0"]][[names(absrel.results[["branch attributes"]][["0"]])[j]]][["Corrected P-value"]]
  }
}

absrel.df <- read.csv("absrel_nodes.csv")
absrel.df <- rbind(absrel.df, data.frame(node=1:length(phy$tip.label), absrel=phy$tip.label))

absrel.df$absrel[which(absrel.df$absrel %in% metadata$name)] <- metadata$tip[match(absrel.df$absrel[which(absrel.df$absrel %in% metadata$name)], metadata$name)]

absrel.p.count <- lengths(lapply(absrel.p, function(x) which(x <=0.05)))
absrel.p.names <- unlist(lapply(absrel.p, function(x) paste(names(which(x <=0.05)), collapse=" ")))
absrel.p.names <- absrel.p.names[absrel.p.names != ""]
absrel.p.names <- gsub("OG000", "", absrel.p.names)

absrel.df$num <- absrel.p.count[match(absrel.df$absrel, names(absrel.p.count))]
absrel.df$orthos <- absrel.p.names[match(absrel.df$absrel, names(absrel.p.names))]
absrel.df$orthos[absrel.df$orthos == ""] <- NA
absrel.df$orthos <- gsub('(?=(?:.{10})+$)', "\n", absrel.df$orthos, perl = TRUE)


gg.selection <- gg.asr %<+% absrel.df +
  geom_tippoint(aes(colour=lifestyle.hyp1), size=5) +
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
  geom_tiplab(fontface=tiplabel.face, size=5, offset=0.2, colour="black") +
  geom_label2(aes(x=branch, label=orthos, fill=num), size=2, label.size=0, fontface="bold") +
  scale_fill_gradient(low="#F0E442", high="#CC79A7",
                      guide=guide_colourbar(title="Effectors showing positive selection",
                                            title.position="top")) +
  theme(legend.box="horizontal",
        legend.title=element_text(size=10, face="bold"),
        legend.text=element_text(size=10))

tiff(file=paste0("fulltree-", Sys.Date(), ".tiff"), height=7, width=12, unit="in", res=300)
gg.selection
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



##MEME

species.tree <- read.tree("../../Phylogenomic analysis/Proteins/RAxML/RAxML_bipartitions.fus_proteins_27T_bsadd")
species.tree <- root(species.tree, "Ilysp1_GeneCatalog_proteins_20121116", edgelabel=TRUE)

meme.gene.orthos <- list()

for (i in rownames(effector.count.SC.SP[effector.count.SC.SP$secretome == "core",])) {
  
  meme.results <- fromJSON(paste0("MEME_results/", i, "_MEME.json"))
  meme.results.gene <- fromJSON(paste0("MEME_results/", i, "_MEME_genetree.json"))
  selection.sites <- which(meme.results$MLE$content$`0`[,7] <= 0.05)
  selection.sites.gene <- which(meme.results.gene$MLE$content$`0`[,7] <= 0.05)
  
  if (all.equal(selection.sites, selection.sites.gene) != TRUE) {
    print(paste(i, length(selection.sites), "species tree |", length(selection.sites.gene), "gene tree"))
  }
  
  gene.tree <- read.tree(paste0("MEME_results/RAxML_bipartitions.", i, "_rooted"))
  gene.tree$tip.label <- metadata$tip[match(gene.tree.5992$tip.label, metadata$file)]
  
  if (length(selection.sites) != 0 & length(selection.sites.gene) != 0) {
    
    meme.gene.orthos[i] <- i
    
    alignment <- fa_read(paste0("MEME_results/", i, "_aln_nuc.translated"))
    
    names(alignment) <- metadata$name[match(names(alignment), metadata$file)]
    alignment <- alignment[tip.order,]
    
    lines <- data.frame(x=0, xend=unique(width(alignment))+1, y=0:length(names(alignment))+0.5, yend=0:length(names(alignment))+0.5)
    rect <- data.frame(xmin=0, xmax=unique(width(alignment))+1, ymin=0.5, ymax=length(names(alignment))+0.5)
    sites <- data.frame(site=selection.sites)
    sites.gene <- data.frame(site=selection.sites.gene)
    
    gg.msa <- ggmsa(alignment, font=NULL, border=FALSE, color="LETTER", posHighligthed=selection.sites) +
      geom_segment(data=lines, aes(x=x, xend=xend, y=y, yend=yend), colour="snow2", size=0.2) +
      geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, colour="dimgrey", size=0.2) +
      geom_text_repel(data=sites, aes(label=site, x=site, y=29), 
                      colour="dimgrey", size=2, angle=90, 
                      #segment.color=NA,
                      force=0.1) +
      ggtitle(sub("OG000", "", i)) +
      scale_x_continuous(breaks=selection.sites,
                         expand=c(0, 0),
                         limits=c(0, 360),
                         position="top") +
      coord_cartesian(clip="off", ylim=c(0, 30)) +
      theme(axis.text.x=element_text(size=3, angle=90, vjust=0.2, hjust=0),
            axis.text.y=element_text(size=4, face=label.face, colour=label.cols),
            plot.title.position="plot",
            plot.title=element_text(size=6, face="bold", vjust=-5),
            plot.margin=margin(0, 0, 5, 0))  
    
    gg.msa.gene <- ggmsa(alignment, font=NULL, border=FALSE, color="LETTER", posHighligthed=selection.sites.gene) +
      geom_segment(data=lines, aes(x=x, xend=xend, y=y, yend=yend), colour="snow2", size=0.2) +
      geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, colour="dimgrey", size=0.2) +
      geom_text_repel(data=sites.gene, aes(label=site, x=site, y=27.5),
                      colour="dimgrey", size=2, 
                      angle=90, 
                      segment.color=NA,
                      direction="x",
                      force=0.1, nudge_y=10) +
      ggtitle(sub("OG000", "", i)) +
      scale_x_continuous(breaks=selection.sites.gene,
                         expand=c(0, 0),
                         limits=c(0, 360),
                         position="top") +
      coord_cartesian(clip="off", ylim=c(0, 32)) +
      theme(axis.text.y=element_text(size=4, face=label.face, colour=label.cols),
            axis.text.x=element_blank(),
            #axis.text.x=element_text(size=3, angle=90, vjust=0.2, hjust=0),
            plot.title=element_text(size=10, face="bold", vjust=-6, hjust=0.1),
            plot.title.position="plot",
            plot.margin=margin(0, 0, 0, 0)) 

    assign(paste0("alignment.", i), alignment)
    assign(paste0("gg.msa.", i), gg.msa)
    assign(paste0("gg.msa.gene.", i), gg.msa.gene)
    
  }
  
  assign(paste0("meme.results.", i), meme.results)
  assign(paste0("meme.results.gene.", i), meme.results.gene)
  assign(paste0("selection.sites.", i), selection.sites)
  assign(paste0("selection.sites.gene", i), selection.sites.gene)
  
}

msa.to.plot <- c("OG0005992", "OG0007067", "OG0006192")


gg.msa.gene <- do.call(ggarrange, c(mget(paste0("gg.msa.gene.", msa.to.plot)), ncol=1))

tiff(file=paste0("test-", Sys.Date(), ".tiff"), height=5, width=10, unit="in", res=300)

gg.msa.gene

dev.off()


do.call(ggarrange, c(mget(paste0("gg.msa.gene.", rownames(effector.count.SC.SP[effector.count.SC.SP$secretome == "core",]))), ncol=1))

#OG with no MEME sites
names(mget(ls(pattern="selection.sites.")))[lengths(mget(ls(pattern="selection.sites."))) == 0]


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





library(ggseqlogo)
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

