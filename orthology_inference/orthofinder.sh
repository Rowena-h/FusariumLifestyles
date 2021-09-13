#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 10    	    # Request 10 cores
#$ -l h_rt=48:00:0      # Request 48 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -m bea

module load anaconda3
conda activate orthofinder

ulimit -n 3944

orthofinder -f . -t 9
