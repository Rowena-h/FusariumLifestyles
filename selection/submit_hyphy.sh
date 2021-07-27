#!/bin/sh

ORTHO=$(cat ../phylogenomics/aln_list | sed 's/\.fa//')

#Number of samples
NUM_ORTHO=$(cat ../phylogenomics/aln_list | wc -l)

mkdir hyphy hyphy/busted hyphy/meme hyphy/absrel hyphy/contrast-fel

#Translate nucleotide alignments for plotting sites under selection
module load emboss

for i in $ORTHO
do
	transeq alignments/codon/${i}_aln_nuc.fa alignments/codon/${i}_aln_nuc.translated
	sed -i 's/_1$//' alignments/codon/${i}_aln_nuc.translated

	#Label gene tree foreground branches for HyPhy
	sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' trees/${i}.raxml.bestTree_rooted > trees/${i}.raxml.bestTree_rooted_hyphy
done

#Label species tree foreground branches for HyPhy
sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' ../divergence_time_estimation/fus_proteins_62T_iqtree_genepart.contree_rooted > trees/fus_proteins_62T_iqtree_genepart.contree_rooted_absrel

module load R

Rscript label_trees.r

qsub -t 1-${NUM_ORTHO} hyphy/BUSTED.sh
qsub -t 1-${NUM_ORTHO} hyphy/aBSREL.sh
qsub -t 1-${NUM_ORTHO} hyphy/Contrast-FEL.sh
