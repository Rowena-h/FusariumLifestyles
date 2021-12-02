#!/bin/sh

STRAINS=$(cat ../strains)

module load abyss

for STRAIN in $STRAINS
do
	abyss-fac abyss/fusotu${STRAIN}/k*/fusotu${STRAIN}-contigs.fa > abyss/fusotu${STRAIN}/fusotu${STRAIN}_abyss_k_comparison.tsv

	KMER=$(tail -n +2 abyss/fusotu${STRAIN}/fusotu${STRAIN}_abyss_k_comparison.tsv | sort -n -k 6 | awk '{print $11}' | sed "s#abyss/fusotu${STRAIN}/##" | sed "s#/fusotu${STRAIN}-contigs.fa##" | tail -n1)

	echo "${KMER} selected" >> abyss/fusotu${STRAIN}/fusotu${STRAIN}_abyss_k_comparison.tsv

	cp abyss/fusotu${STRAIN}/${KMER}/fusotu${STRAIN}-contigs.fa abyss/fusotu${STRAIN}/
done
