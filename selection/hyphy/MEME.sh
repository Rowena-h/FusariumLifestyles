#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 cores
#$ -l h_rt=1:00:00 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -t 18
#$ -j y

ORTHO=$(cat ../orthogroups_selection.csv | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate hyphy

hyphy meme --alignment ../alignments/${ORTHO}_aln_nuc.fa --tree RAxML_bipartitions.fus_proteins_27T_bsadd_rooted --output raxml/${ORTHO}_MEME.json
