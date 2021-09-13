#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1 		# Request 1 core
#$ -l h_rt=1:00:00 	# Request 1 hour runtime
#$ -l h_vmem=3G   	# Request 3GB RAM
#$ -j y
#$ -m bea
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load samtools

samtools index ../polishing/fusotu${STRAIN}_abyss_mapped_sorted_dups.bam

module load anaconda3
conda activate blobtools

~/Programmes/blobtools/blobtools create 	-i ../denovo_assembly/abyss/fusotu${STRAIN}/fusotu${STRAIN}-contigs.fa \
						-b ../polishing/fusotu${STRAIN}_abyss_mapped_sorted_dups.bam \
						-t fusotu${STRAIN}_abyss_blast.tsv \
						-o fusotu${STRAIN}_abyss

~/Programmes/blobtools/blobtools plot 	-r species \
					-i fusotu${STRAIN}_abyss.blobDB.json
