#!/bin/sh
#$ -cwd
#$ -pe parallel 48
#$ -l infiniband=sdv-i
#$ -l h_rt=24:0:0
#$ -m bea
#$ -t 120-128	#kmer sizes to check

module load abyss

for STRAIN in 1 3 5 6 7
do
	mkdir abyss/fusotu${STRAIN}/k${SGE_TASK_ID}
	abyss-pe -C abyss/fusotu${STRAIN}/k${SGE_TASK_ID} name=fusotu${STRAIN} in='FUS_OTU${STRAIN}_1_trimmedpaired.fastq.gz FUS_OTU${STRAIN}_2_trimmedpaired.fastq.gz' np=$NSLOTS
done
