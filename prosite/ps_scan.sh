#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1    	# Request 48 cores
#$ -l h_rt=01:00:00     # Request 72 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -m bea

DIR=$1
SAMPLES=$(cat $2)
RUN_DIR=/data/SBCS-BuggsLab/RowenaHill/fus_comparison

for i in $SAMPLES
do
	perl /data/home/btx494/Programmes/ps_scan/ps_scan.pl -p PS00014 -o scan -d /data/home/btx494/Programmes/ps_scan/prosite.dat ${DIR}${i} > ${RUN_DIR}/prosite/${i}_psscan
	grep ">" ${RUN_DIR}/prosite/${i}_psscan | awk '{print $1}' | tr -d '>' > ${RUN_DIR}/prosite/${i}_psscan_ERlist
done
