#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 core
#$ -l h_rt=12:00:00 	# Request 30 minutes runtime
#$ -l h_vmem=20G   	# Request 1GB RAM
#$ -j y

TAXON=$(ls ../../proteins/*.faa | sed 's#\.\./\.\./proteins/##' | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate run_dbcan

run_dbcan ../../proteins/${TAXON} protein --out_dir ${TAXON}_dbcan_results

awk '$6 == 3' ${TAXON}_dbcan_results/overview.txt > ${TAXON}_cazymes
