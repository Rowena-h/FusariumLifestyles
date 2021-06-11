#!/bin/sh

ORTHO=$(ls -1 ../../phylogenomics/gene_trees/*_aln.fa | sed 's#\.\./\.\./phylogenomics/gene_trees/##' | sed 's/_aln\.fa//')
CORE_EFFECTORS=$(cat orthogroups_meme.csv)
#Number of samples
NUM_ORTHO=$(ls -1 ../../phylogenomics/gene_trees/*_aln.fa | wc -l)
NUM_CORE_EFFECTORS=$(cat orthogroups_meme.csv | wc -l)

#mkdir hyphy hyphy/busted hyphy/meme hyphy/absrel codon

#Translate nucleotide alignments for plotting sites under selection
#module load emboss

#for i in $ORTHO
#do
#	transeq ls ../divergence_time_estimation/alignments/codon/${i}_aln_nuc.fa codon/${i}_aln_nuc.translated
#	sed -i 's/_1$//' codon/${i}_aln_nuc.translated

#	#Label gene tree foreground branches for HyPhy
#	sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' ../divergence_time_estimation/trees/RAxML_bipartitions.${ORTHO}_rooted
#done

#Label species tree foreground branches for HyPhy
#sed -r 's/:/\{FOREGROUND\}:/g; s/\{FOREGROUND\}([^,]*)$/\1/' ../divergence_time_estimation/fus_proteins_62T_iqtree_genepart.contree_rooted fus_proteins_62T_iqtree_genepart.contree_rooted_absrel

#qsub -t 1-${NUM_ORTHO} BUSTED.sh
#qsub -t 1-${NUM_CORE_EFFECTORS} MEME.sh
qsub -t 1-${NUM_ORTHO} aBSREL.sh
