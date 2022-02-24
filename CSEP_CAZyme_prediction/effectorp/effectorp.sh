#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1    	    # Request 1 core
#$ -l h_rt=00:30:0      # Request 30 minutes runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

DIR=$1
SAMPLE=$(cat $2 | sed -n ${SGE_TASK_ID}p)
RUN_DIR=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/CSEP_prediction

python /data/home/btx494/Programmes/EffectorP-3.0/EffectorP.py -i ${DIR}${SAMPLE} > ${RUN_DIR}/effectorp/${SAMPLE}_effectorp
#List of effectors
awk '/# Identifier/{flag=1;next}/-----------------/{flag=0}flag' ${RUN_DIR}/effectorp/${SAMPLE}_effectorp | awk '$NF=="effector" { print $1}' | sed '/^$/d' > ${RUN_DIR}/effectorp/${SAMPLE}_effectorlist
