#!/bin/sh

ORTHO=$(ls -1 ../phylogenomics/gene_trees/*_aln.fa | sed 's#\.\./phylogenomics/gene_trees/##' | sed 's/_aln\.fa//')
CORE_EFFECTORS=$(cat ../orthogroups_selection.csv)
#Number of samples
NUM_ORTHO=$(ls -1 ../phylogenomics/gene_trees/*_aln.fa | sed 's#\.\./phylogenomics/gene_trees/##' | sed 's/_aln\.fa//' | wc -l)
NUM_CORE_EFFECTORS=$(cat ../orthogroups_selection.csv | wc -l)

module load emboss

for i in $ORTHO
do
	transeq ../alignments/${i}_aln_nuc.fa ../alignments/${i}_aln_nuc.translated
	sed -i 's/_1$//' ../alignments/${i}_aln_nuc.translated
done

qsub -t 1-${NUM_ORTHO} BUSTED.sh
#qsub -t 1-${NUM_CORE_EFFECTORS} aBSREL.sh MEME.sh MEME_gene.sh
