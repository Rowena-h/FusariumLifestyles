#!/bin/sh
#Script to submit de novo assembly jobs

for ASSEMBLER in abyss megahit spades
do
	mkdir ${ASSEMBLER}

	for i in 1 3 5 6 7
	do
		mkdir ${ASSEMBLER}/fusotu${i}
	done

	qsub ${ASSEMBLER}.sh
done
