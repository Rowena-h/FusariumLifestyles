#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1        	# Request 1 core
#$ -l h_rt=1:00:0 	# Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -t 2-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../../strains)
TAG=$(sed -n ${SGE_TASK_ID}p locus_tags.txt)

cd fusotu${STRAIN}_abyss_rnd3.maker.output

module load maker

maker_map_ids 	--prefix ${TAG}_ \
		--justify 6 fusotu${STRAIN}_abyss_rnd3.all.maker.gff > fusotu${STRAIN}_abyss_rnd3.all.maker.map

map_gff_ids fusotu${STRAIN}_abyss_rnd3.all.maker.map fusotu${STRAIN}_abyss_rnd3.all.maker.gff

map_fasta_ids fusotu${STRAIN}_abyss_rnd3.all.maker.map fusotu${STRAIN}_abyss_rnd3.all.maker.proteins.fasta

map_fasta_ids fusotu${STRAIN}_abyss_rnd3.all.maker.map fusotu${STRAIN}_abyss_rnd3.all.maker.transcripts.fasta
