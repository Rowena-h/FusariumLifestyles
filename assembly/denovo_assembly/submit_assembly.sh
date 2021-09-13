#!/bin/sh
#Script to submit de novo assembly jobs

STRAINS=$(cat ../strains)

for ASSEMBLER in abyss megahit spades
do
	mkdir ${ASSEMBLER}

	for STRAIN in $STRAINS
	do
		mkdir ${ASSEMBLER}/fusotu${STRAIN}
	done

	qsub ${ASSEMBLER}.sh
done
