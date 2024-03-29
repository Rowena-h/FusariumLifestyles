#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 core
#$ -l h_rt=5:00:00      # Request 5 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -m bea
#$ -j y

mkdir run${SGE_TASK_ID}_correlated
cp {in.BV,fus_proteins_62T.raxml.support_dating,fus_proteins_dating10_mcmctree_short.phy,mcmctree_step2_correlated.ctl} run${SGE_TASK_ID}_correlated

cd run${SGE_TASK_ID}_correlated
sed -i "s/mcmcfile = mcmc.txt/mcmcfile = mcmc_run${SGE_TASK_ID}_correlated.txt/" mcmctree_step2_correlated.ctl

module load anaconda3
conda activate paml

mcmctree mcmctree_step2_correlated.ctl
