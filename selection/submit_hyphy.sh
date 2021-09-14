#!/bin/sh

ORTHO=$(cat ../phylogenomics/aln_list | sed 's/\.fa//')

#Number of samples
NUM_ORTHO=$(cat ../phylogenomics/aln_list | wc -l)

mkdir hyphy hyphy/busted hyphy/meme hyphy/absrel hyphy/contrast-fel

#Translate nucleotide alignments for plotting sites under selection
#module load emboss

#for i in $ORTHO
#do
#	transeq alignments/codon/${i}_aln_nuc.fa alignments/codon/${i}_aln_nuc.translated
#	sed -i 's/_1$//' alignments/codon/${i}_aln_nuc.translated
#
#	#Label gene tree foreground branches for HyPhy
#	sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' trees/${i}.raxml.bestTree_rooted > trees/${i}.raxml.bestTree_rooted_hyphy
#done

#Label dated species tree foreground branches for aBSREL and BUSTED
sed -n 638p ../divergence_time_estimation/mcmctree/run1_independent/mcmctree_step2_out.txt | sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' > trees/dated_tree_absrel.tre
awk -F',' '{print $8 "\t" $9}' ../metadata.csv > tmp
sed `cat tmp | awk '{print "-e s/"$2"/"$1"/"}'`<<<"`cat trees/dated_tree_absrel.tre`" > tmp2 && mv tmp2 trees/dated_tree_absrel.tre
rm tmp

#Label trees for different lifestyles for Contrast-FEL
module load R

Rscript label_trees.r

#Submit HyPhy programmes
qsub -t 1-${NUM_ORTHO} hyphy/BUSTED.sh
qsub -t 1-${NUM_ORTHO} hyphy/aBSREL.sh
qsub -t 1-${NUM_ORTHO} hyphy/Contrast-FEL.sh
