#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 4 cores
#$ -l h_rt=0:30:0 	# Request 30 min runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -j y

STRAINS=$(cat ../strains)

module load anaconda3
conda activate quast

for STRAIN in $STRAINS
do
	quast.py ../polishing/fusotu${STRAIN}_abyss_pilon_filtered.fa ../polishing/fusotu${STRAIN}_megahit_pilon_filtered.fa ../polishing/fusotu${STRAIN}_spades_pilon_filtered.fa -o fusotu${STRAIN}_quast_results -t ${NSLOTS} --fungus -l "ABySS v2.0.2, MEGAHIT v1.2.9, SPAdes v3.11.1"
done
