#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 cores
#$ -l h_rt=12:00:00 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

ORTHO=$(cat ../../phylogenomics/aln_list | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate hyphy-2.5.30

hyphy busted 	--alignment ..//alignments/codon/${ORTHO}_aln_nuc.fa \
		--tree ../trees/${ORTHO}.raxml.bestTree_rooted \
		--output busted/${ORTHO}_BUSTED.json \
		--branches FOREGROUND
