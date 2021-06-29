#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 5            # Request 1 cores
#$ -l h_rt=1:00:00      # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

ORTHO=$(cat aln_list2 | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)

MODEL=$(grep ${ORTHO} species_tree/fus_proteins_62T_raxmlpartition.txt | sed 's/,.*$//')

module load anaconda3
conda activate raxml-ng

raxml-ng --search \
	 --msa gene_trees/${ORTHO}_aln_edit_trimmed.phy \
	 --model ${MODEL} \
	 --prefix gene_trees/RAxML-NG/${ORTHO} \
	 --seed 2 \
	 --threads ${NSLOTS}

#raxml-ng --bsconverge \
#        --bs-trees gene_trees/RAxML-NG/${ORTHO}.raxml.bootstraps \
#        --prefix gene_trees/RAxML-NG/${ORTHO}_convergence_test \
#        --seed 2
