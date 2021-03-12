#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1    	# Request 48 cores
#$ -l h_rt=24:00:0     # Request 72 hours runtime
#$ -l h_vmem=3G         # Request 1GB RAM
#$ -j y

DIR=$1
SAMPLE=$(cat $2 | sed -n ${SGE_TASK_ID}p)
RUN_DIR=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/effector_prediction

/data/home/btx494/Programmes/targetp-2.0/bin/targetp -fasta ${DIR}${SAMPLE} -org non-pl -prefix ${RUN_DIR}/targetp/${SAMPLE}_targetp -tmp ${RUN_DIR}/targetp/tmp_${SAMPLE}
#List of signal peptide genes
cat ${RUN_DIR}/targetp/${SAMPLE}_targetp_summary.targetp2 | awk '$2=="SP" { print $1}' > ${RUN_DIR}/targetp/${SAMPLE}_targetp_SPlist
