#!/bin/sh
#$ -cwd  		    # Set the working directory for the job to the current directory
#$ -pe smp 1    	# Request 1 core
#$ -l h_rt=24:00:00 # Request 24 hours runtime
#$ -l h_vmem=5G     # Request 5GB RAM
#$ -j y

DIR=$1
SAMPLES=$(cat $2)

cd /data/home/btx494/Programmes/nucpred-1.1/

for i in $SAMPLES
do
	/data/home/btx494/Programmes/nucpred-1.1/nucpred-rh.pl ${DIR}${i} > ${i}_nucpred
	cat ${i}_nucpred | awk '$NF>=0.8 {print $1}' > ${i}_nucpred_list

        PROTEINS=$(grep ">" /data/SBCS-BuggsLab/RowenaHill/fus_comparison/orthology_inference/${i} | wc -l)
        LENGTH=$(cat ${i}_nucpred | wc -l)

        if [ "$PROTEINS" -ne "$LENGTH" ]
        then
                echo $i >> nucpred_failed
        else
		mv ${i}_nucpred* /data/SBCS-BuggsLab/RowenaHill/fus_comparison/CSEP_prediction/nucpred/
	fi
done

mv nucpred_failed /data/SBCS-BuggsLab/RowenaHill/fus_comparison/CSEP_prediction/nucpred/
