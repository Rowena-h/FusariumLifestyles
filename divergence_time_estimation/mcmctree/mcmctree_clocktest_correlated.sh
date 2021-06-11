#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 cores
#$ -l h_rt=240:00:00     # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -m bea
#$ -j y

cd clock_testing/correlated/
sed -i "s/clock = /clock = 3/" mcmctree_clocktest_correlated.ctl

module load anaconda3
conda activate paml

mcmctree mcmctree_clocktest_correlated.ctl
