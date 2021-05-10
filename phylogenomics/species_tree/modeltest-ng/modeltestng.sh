#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 6            # Request 1 cores
#$ -l h_rt=240:00:00      # Request 1 hour runtime
#$ -l h_vmem=10G         # Request 1GB RAM
#$ -j y
#$ -m bea

PARTITION=$(ls -dir partition.* | sed -n ${SGE_TASK_ID}p)

cd partition.${SGE_TASK_ID}

module load anaconda3
conda activate modeltest-ng

modeltest-ng -d aa -i fus_proteins_62T_concat.phy -q fus_proteins_62T_partition.num${SGE_TASK_ID} -p ${NSLOTS} -T raxml
