#!/bin/sh

#Make fasta file of core single copy orthogroups for each taxon

for ORTHO in $(ls ../alignments/codon/*.fa | sed 's#\.\./alignments/codon/##')
do
	awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' ../alignments/codon/${ORTHO} | sed 's#\.\./alignments/codon/##' > ${ORTHO}.tmp
done

for TAXON in $(ls ../../proteins/*.faa | sed 's#\.\./\.\./proteins/##' | sed 's/\.faa//')
do
	awk -v p="$TAXON" 'BEGIN{ ORS=""; RS=">"; FS="\n" } $1 ~ p { print ">" $0 }' *.tmp | sed "s/_${TAXON}//" | sed 's/_aln_nuc\.fa//' > ${TAXON}_coreSC.fa
done

rm *.tmp

module load R/4.0.2

Rscript codon_optimisation.r ../../orthology_inference/OrthoFinder/Results_Oct22/
