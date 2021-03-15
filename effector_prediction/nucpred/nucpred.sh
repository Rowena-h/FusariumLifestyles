#!/bin/sh
#$ -cwd  		# Set the working directory for the job to the current directory
#$ -pe smp 1    	# Request 48 cores
#$ -l h_rt=12:00:00     # Request 72 hours runtime
#$ -l h_vmem=2G         # Request 1GB RAM
#$ -j y
#$ -m bea

DIR=$1
SAMPLES=$(cat $2)

cd /data/home/btx494/Programmes/nucpred-1.1/

for i in $SAMPLES
do
	/data/home/btx494/Programmes/nucpred-1.1/nucpred-rh.pl ${DIR}${i} > ${i}_nucpred
	cat ${i}_nucpred | awk '$NF>=0.8 {print $1}' > ${i}_nucpred_list
	#mv ${i}_nucpred* /data/SBCS-BuggsLab/RowenaHill/fus_comparison/effector_prediction/nucpred/
done
