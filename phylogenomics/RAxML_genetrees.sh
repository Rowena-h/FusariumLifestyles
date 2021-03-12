#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 cores
#$ -l h_rt=72:00:00 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y
#$ -t 1

ORTHO=$(cat aln_list | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)
ORTHOFINDER_DIR=../orthology_inference/OrthoFinder/Results_Mar01/Single_Copy_Orthologue_Sequences

#FIX NAMES

sed -e 's/.faa.*//' ${ORTHOFINDER_DIR}/${ORTHO}.fa > gene_trees/${ORTHO}_edit.fa

#ALIGN

module load mafft

mafft gene_trees/${ORTHO}_edit.fa > gene_trees/${ORTHO}_aln.fa

#TRIM

java -jar /data/home/btx494/Programmes/BMGE-1.12/BMGE.jar -i gene_trees/${ORTHO}_aln.fa -t AA -o gene_trees/${ORTHO}_alntrimmed.phy 

#MAKE GENE TREES

module load raxml

raxmlHPC-SSE3 	-f a \
		-n ${ORTHO} \
		-s gene_trees/${ORTHO}_alntrimmed.phy \
		-N 1000 \
		-m PROTGAMMAAUTO \
		-p 12345 \
		-x 12345 \
		-w /data/SBCS-BuggsLab/RowenaHill/fus_comparison/phylogenomics/gene_trees/RAxML

rm gene_trees/${ORTHO}_alntrimmed.phy.reduced

