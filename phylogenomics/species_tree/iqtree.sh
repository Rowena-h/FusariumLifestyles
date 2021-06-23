#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 8            # Request 1 cores
#$ -l h_rt=240:00:00      # Request 1 hour runtime
#$ -l h_vmem=5G         # Request 1GB RAM
#$ -j y
#$ -m bea

module load anaconda3
conda activate IQ-Tree

iqtree -s fus_proteins_62T_concat.phy -spp fus_proteins_62T_iqtreepartition.txt -bb 1000 -nt ${NSLOTS} -pre iqtree/fus_proteins_62T_iqtree_genepart -m MFP
