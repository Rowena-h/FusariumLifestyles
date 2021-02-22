#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 cores
#$ -l h_rt=24:00:00 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y
#$ -o /dev/null
#$ -t 1-1064

ORTHO=$(cat aln_list | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)

#FIX NAMES

sed -e 's/.faa.*//' ${ORTHO}.fa > ${ORTHO}_edit.fa

#ALIGN

module load mafft

mafft ${ORTHO}_edit.fa > ${ORTHO}_aln.fa

#TRIM

java -jar /data/home/btx494/Programmes/BMGE-1.12/BMGE.jar -i ${ORTHO}_aln.fa -t DNA -o ${ORTHO}_alntrimmed.phy 

#MAKE GENE TREES

module load raxml

mkdir trees

raxmlHPC-SSE3 	-f a \
		-n ${ORTHO} \
		-s ${ORTHO}_alntrimmed.phy \
		-N 1000 \
		-m GTRGAMMA \
		-p 12345 \
		-x 12345 \
		-w /data/SBCS-BuggsLab/RowenaHill/fus_comparison/orthofinder/OrthoFinder/Results_Feb19/Single_Copy_Orthologue_Sequences/trees

rm ${ORTHO}_alntrimmed.phy.reduced

