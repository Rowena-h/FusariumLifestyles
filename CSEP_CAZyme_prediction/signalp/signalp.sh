#!/bin/sh
#$ -cwd 		    # Set the working directory for the job to the current directory
#$ -pe smp 1    	# Request 1 core
#$ -l h_rt=48:00:0  # Request 48 hours runtime
#$ -l h_vmem=5G     # Request 5GB RAM
#$ -j y
#$ -m bea

DIR=$1
SAMPLES=$(cat $2)

cd /data/home/btx494/Programmes/signalp-5.0b/bin/

for i in $SAMPLES
do
	./signalp -fasta ${DIR}${i} -org euk -prefix ${i}_signalp
	#List of signal peptide genes
        cat ${i}_signalp_summary.signalp5 | awk '$2=="SP(Sec/SPI)" {print $1}' > ${i}_signalp_SPlist
	mv ${i}_signalp* /data/SBCS-BuggsLab/RowenaHill/fus_comparison/CSEP_prediction/signalp/
done
