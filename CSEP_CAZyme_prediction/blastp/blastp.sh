#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 core
#$ -l h_rt=1:00:00 	# Request 30 minutes runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

TAXON=$(ls ../../proteins/*.faa | sed 's#\.\./\.\./proteins/##' | sed 's/\.faa//' | sed -n ${SGE_TASK_ID}p)

module load seqtk

seqtk subseq ../../proteins/${TAXON}.faa ../${TAXON}.faa_candidate_effectors > ${TAXON}_cseps.faa

module load blast+

blastp \
-query ${TAXON}_cseps.faa \
-db phi-base_current.fas \
-outfmt '6 qseqid sseqid evalue bitscore pident length' \
-evalue 1e-25 \
-out ${TAXON}_phi-b_blast \
-num_threads ${NSLOTS}

#Filter for top bitscore result per gene
sort -r -n -k4 < ${TAXON}_phi-b_blast | awk '!x[$1]++' ${TAXON}_phi-b_blast | sort -k 1b,1  > ${TAXON}_phi-b_blast_tophits
#Add data to CSEPs file
join -a1 -a2 -t $'\t' -o 1.1 2.2 2.3 2.4 2.5 -1 1 -2 1 ../${TAXON}.faa_candidate_effectors ${TAXON}_phi-b_blast_tophits > ../${TAXON}.faa_cseps
awk -F "\t" '{gsub(/\#/,"\t",$2);print $0}' ../${TAXON}.faa_cseps | sed 's/ /\t/g' > ${TAXON}tmp && mv ${TAXON}tmp ../${TAXON}.faa_cseps
