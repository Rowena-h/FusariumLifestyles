#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 cores
#$ -l h_rt=00:30:00     # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

module load anaconda3
conda activate phyx

python get_var_length.py ../../../selection/trees/ --flend _rooted --outf var --outg Ilysp1_GeneCatalog_proteins_20121116

python get_bp_genetrees.py ../../../selection/trees/ ../../../selection/hyphy/absrel/fus_proteins_62T_iqtree_genepart.contree_rooted --flend _rooted --outf bp

python combine_results.py var bp --outf comb

python get_good_genes.py comb --max 10 --outf dating_orthogroups
