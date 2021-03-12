#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1    	# Request 48 cores
#$ -l h_rt=12:00:00     # Request 72 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -m bea

DIR=$1
SAMPLES=$(cat $2)
RUN_DIR=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/effector_prediction

for i in $SAMPLES
do
	/data/home/btx494/Programmes/tmhmm-2.0c/bin/tmhmm ${DIR}${i} > ${RUN_DIR}/tmhmm/${i}_tmhmm
	grep "Number of predicted TMHs" ${RUN_DIR}/tmhmm/${i}_tmhmm | awk '{ $NF = "\t" $NF; print }' | column -t -s $'\t' | awk '$7>1 { print $2}' > ${RUN_DIR}/tmhmm/${i}_tmhmm_TMlist
done
