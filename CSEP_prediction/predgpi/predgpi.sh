#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1    	    # Request 1 core
#$ -l h_rt=24:00:0      # Request 24 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -m bea

DIR=$1
SAMPLES=$(cat $2)

module load R

Rscript PredGPI.r ${DIR} ${SAMPLES}