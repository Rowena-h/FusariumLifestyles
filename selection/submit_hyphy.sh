#!/bin/sh

#ORTHO=$(cat ../phylogenomics/aln_list | sed 's/\.fa//')
#CORE_EFFECTORS=$(cat orthogroups_meme.csv)
#ORTHO=$(cat hyphy/torepeat2)

#Number of samples
#NUM_ORTHO=$(cat ../phylogenomics/aln_list | wc -l)
#NUM_CORE_EFFECTORS=$(cat orthogroups_meme.csv | wc -l)
#NUM_ORTHO=$(cat hyphy/torepeat2 | wc -l)

#mkdir trees hyphy/busted hyphy/meme hyphy/absrel hyphy/contrast-fel codon

#Translate nucleotide alignments for plotting sites under selection
#module load emboss

#for i in $ORTHO
#do
	#transeq ../divergence_time_estimation/alignments/codon/${i}_aln_nuc.fa codon/${i}_aln_nuc.translated
	#sed -i 's/_1$//' codon/${i}_aln_nuc.translated

	#Label gene tree foreground branches for HyPhy
	#sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' ../divergence_time_estimation/trees/${i}.raxml.bestTree_rooted > trees/${i}.raxml.bestTree_rooted
#done

#Label species tree foreground branches for HyPhy
#sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' ../divergence_time_estimation/fus_proteins_62T_iqtree_genepart.contree_rooted > trees/fus_proteins_62T_iqtree_genepart.contree_rooted_absrel

module load R

Rscript label_trees.r

#qsub -t 1-${NUM_ORTHO} BUSTED.sh
#qsub -t 1-${NUM_CORE_EFFECTORS} MEME.sh
#qsub -t 1-${NUM_ORTHO} aBSREL.sh
#qsub -t 1-${NUM_ORTHO} Contrast-FEL.sh
