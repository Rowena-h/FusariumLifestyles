#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 5            # Request 5 cores
#$ -l h_rt=1:00:00      # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

ORTHO=$(cat aln_list | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)

MODEL=$(grep ${ORTHO} species_tree/fus_proteins_62T_raxmlpartition.txt | sed 's/,.*$//')
MODEL_BMGE=$(grep ${ORTHO} species_tree/fus_proteins_bmge_62T_raxmlpartition.txt | sed 's/,.*$//')

module load anaconda3
conda activate raxml-ng

raxml-ng --search \
         --msa gene_trees/${ORTHO}_aln_edit_trimmed.phy \
         --model ${MODEL} \
         --prefix gene_trees/RAxML-NG/${ORTHO} \
         --seed 2 \
         --threads ${NSLOTS}

raxml-ng --search \
         --msa gene_trees_bmge/${ORTHO}_aln_edit_trimmed.phy \
         --model ${MODEL_BMGE} \
         --prefix gene_trees_bmge/RAxML-NG/${ORTHO} \
         --seed 2 \
         --threads ${NSLOTS}
