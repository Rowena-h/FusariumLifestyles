#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=48:00:00 # Request 48 hours runtime
#$ -l h_vmem=2G   	# Request 2GB RAM
#$ -m bea
#$ -j y
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load busco

export AUGUSTUS_CONFIG_PATH=/data/SBCS-BuggsLab/RowenaHill/genome_assemblies/augustus_config/config

for ASSEMBLER in abyss megahit spades
do
	BUSCO.py -i ../polishing/fusotu${STRAIN}_${ASSEMBLER}_pilon.fasta -c ${NSLOTS} -o fusotu${STRAIN}_${ASSEMBLER}_pilon -m genome -l /data/SBCS-BuggsLab/RowenaHill/genome_assemblies/busco_datasets/hypocreales_odb10.2019-11-20 -sp fusarium
done
