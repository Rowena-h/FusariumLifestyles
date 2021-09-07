#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 20 cores
#$ -l h_rt=24:00:0 	# Request 24 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load anaconda3
conda activate repeatmasker

mkdir fusotu${STRAIN}abyss_masked
RepeatMasker -e ncbi -lib fusotu${STRAIN}/RM*/consensi.fa -pa ${NSLOTS} -xsmall -dir fusotu${STRAIN}_abyss_masked ../../assembly/polishing/fusotu${STRAIN}_abyss_pilon.fasta
