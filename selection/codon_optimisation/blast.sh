#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1 		# Request 4 cors
#$ -l h_rt=1:00:00 	# Request 30 min runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y 

TAXON=$(ls ../../proteins/*.faa | sed 's#\.\./\.\./proteins/##' | sed -n ${SGE_TASK_ID}p)

module load blast+

makeblastdb -dbtype prot -in ../../proteins/${TAXON}

blastp \
-query ribosomal_proteins.fasta \
-db ../../proteins/${TAXON} \
-outfmt '6 qseqid sseqid evalue bitscore pident length' \
-evalue 1e-25 \
-out ${TAXON}_blast \
-num_threads ${NSLOTS}

#Make file with list of protein names with blast hits
awk '{print $2}' ${TAXON}_blast > ${TAXON}_ribosomes
#grep -f ${TAXON}_blast_ribosomes ../../orthology_inference/OrthoFinder/Results_Mar01/Orthogroups/Orthogroups.tsv | awk '{print $1}' > ${TAXON}_ribosomeorthos
