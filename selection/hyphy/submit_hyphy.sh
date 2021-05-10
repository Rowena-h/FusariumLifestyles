#!/bin/sh

ORTHO=$(ls -1 ../../phylogenomics/gene_trees/*_aln.fa | sed 's#\.\./\.\./phylogenomics/gene_trees/##' | sed 's/_aln\.fa//')
CORE_EFFECTORS=$(cat ../orthogroups_selection.csv)
#Number of samples
NUM_ORTHO=$(ls -1 ../../phylogenomics/gene_trees/*_aln.fa | wc -l)
NUM_CORE_EFFECTORS=$(cat ../orthogroups_selection.csv | wc -l)

mkdir busted meme absrel

#Translate nucleotide alignments for plotting sites under selection
module load emboss

for i in $CORE_EFFECTORS
do
	transeq ../alignments/codon/${i}_aln_nuc.fa ../alignments/codon/${i}_aln_nuc.translated
	sed -i 's/_1$//' ../alignments/codon/${i}_aln_nuc.translated
done

qsub -t 1-${NUM_ORTHO} BUSTED.sh
qsub -t 1-${NUM_CORE_EFFECTORS} MEME_gene.sh
qsub -t 1-${NUM_CORE_EFFECTORS} aBSREL.sh
