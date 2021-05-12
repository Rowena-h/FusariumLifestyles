#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1           # Request 48 cores
#$ -l h_rt=00:30:00     # Request 72 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

module load anaconda3
conda activate endophyte_determinants

python svm.py
