#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 cores
#$ -l h_rt=5:00:00 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

ORTHO=$(cat ../../phylogenomics/aln_list | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate hyphy-2.5.30

hyphy absrel 	--alignment ../alignments/codon/${ORTHO}_aln_nuc.fa \
		--tree ../trees/fus_proteins_62T_iqtree_genepart.contree_rooted \
		--output absrel/${ORTHO}_aBSREL.json \
		--branches FOREGROUND