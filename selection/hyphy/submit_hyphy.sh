#!/bin/sh

ORTHO=$(cat ../orthogroups_selection.csv)

module load emboss

for i in $ORTHO
do
	transeq ../alignments/${i}_aln_nuc.fa ../alignments/${i}_aln_nuc.translated
	sed -i 's/_1$//' ../alignments/${i}_aln_nuc.translated
done

qsub aBSREL.sh MEME.sh MEME_gene.sh
