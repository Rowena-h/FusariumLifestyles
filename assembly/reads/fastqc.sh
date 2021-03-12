#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 4      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 1GB RAM
#$ -j y
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load fastqc

fastqc -t ${NSLOTS} ${STRAIN}_1_trimmedpaired.fastq.gz ${STRAIN}_2_trimmedpaired.fastq.gz
