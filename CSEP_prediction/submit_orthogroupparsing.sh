#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 48 cores
#$ -l h_rt=1:00:00     # Request 72 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -m bea

module load R

./orthogroup_parser.r /data/SBCS-BuggsLab/RowenaHill/fus_comparison/orthology_inference/OrthoFinder/Results_Mar01/
