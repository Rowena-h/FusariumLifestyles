#!/bin/sh
#$ -cwd
#$ -pe parallel 48
#$ -l infiniband=sdv-i
#$ -l h_rt=24:0:0
#$ -m bea
#$ -t 120-128	#kmer sizes to check

STRAINS=$(cat ../strains)

module load abyss

for STRAIN in $STRAINS
do
	mkdir abyss/fusotu${STRAIN}/k${SGE_TASK_ID}
	abyss-pe -C abyss/fusotu${STRAIN}/k${SGE_TASK_ID} name=fusotu${STRAIN} in='../reads/FUS_OTU${STRAIN}_1_trimmedpaired.fastq.gz ../reads/FUS_OTU${STRAIN}_2_trimmedpaired.fastq.gz' np=${NSLOTS}
done
