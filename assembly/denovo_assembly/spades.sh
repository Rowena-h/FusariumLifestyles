#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 1 core
#$ -l h_rt=48:0:0 	# Request 24 hour runtime
#$ -l highmem
#$ -l h_vmem=60G   	# Request 1GB RAM
#$ -j y
#$ -m bea
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load spades

spades.py 	-1 ../reads/FUS_OTU${STRAIN}_1_trimmedpaired.fastq.gz \
		-2 ../reads/FUS_OTU${STRAIN}_2_trimmedpaired.fastq.gz \
		--careful \
		-o spades/fusotu${STRAIN} \
		-t ${NSLOTS}

mv spades/fusotu${STRAIN}/contigs.fasta spades/fusotu${STRAIN}/fusotu${STRAIN}-contigs.fa
