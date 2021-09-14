#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 core
#$ -l h_rt=1:00:00      # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -m bea

module load R

Rscript orthogroup_parser.r ../orthology_inference/OrthoFinder/Results_Mar01/
