#!/bin/sh

ORTHO=$(cat ../phylogenomics/aln_list | sed 's/\.fa//')

#Number of samples
NUM_ORTHO=$(cat ../phylogenomics/aln_list | wc -l)

mkdir hyphy/busted hyphy/absrel hyphy/contrast-fel

for i in $ORTHO
do
	#Label gene tree foreground branches for HyPhy
	sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' trees/${i}.raxml.bestTree_rooted > trees/${i}.raxml.bestTree_rooted_hyphy
done

#Label dated species tree foreground branches for aBSREL and BUSTED
sed -n 631p ../divergence_time_estimation/mcmctree/run1_independent/mcmctree_step2_out.txt | sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' > trees/dated_tree_absrel.tre
awk -F',' '{print $10 "\t" $11}' ../metadata.csv > tmp
sed `cat tmp | awk '{print "-e s/"$2"/"$1"/"}'`<<<"`cat trees/dated_tree_absrel.tre`" > tmp2 && mv tmp2 trees/dated_tree_absrel.tre
rm tmp

#Label trees for different lifestyles for Contrast-FEL
module load R

Rscript label_trees.r

cd hyphy 

#Submit HyPhy programmes
qsub -t 1-${NUM_ORTHO} BUSTED.sh
qsub -t 1-${NUM_ORTHO} aBSREL.sh
qsub -t 1-${NUM_ORTHO} Contrast-FEL.sh
