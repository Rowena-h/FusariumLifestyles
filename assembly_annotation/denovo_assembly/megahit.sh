#!/bin/sh
#$ -cwd
#$ -pe smp 10
#$ -l h_rt=48:0:0
#$ -m bea
#$ -t 1-5

STRAIN=$(cat strains | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate megahit

megahit -1 FUS_OTU${STRAIN}_1_trimmedpaired.fastq.gz -2 FUS_OTU${STRAIN}_1_trimmedpaired.fastq.gz -t ${NSLOTS} -o megahit/fusotu${STRAIN} --k-list 51,59,67,75,83,91,99,107,115,123,131
