#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 cores
#$ -l h_rt=72:00:00 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

#ORTHO=$(cat aln_list | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)
ORTHO=$(cat aln_list2 | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)
ORTHOFINDER_DIR=../orthology_inference/OrthoFinder/Results_Mar01/Single_Copy_Orthologue_Sequences

#ALIGN

#module load mafft

#mafft ${ORTHOFINDER_DIR}/${ORTHO}.fa > gene_trees/${ORTHO}_aln.fa

#FIX NAMES

#sed 's/.faa.*$//' gene_trees/${ORTHO}_aln.fa > gene_trees/${ORTHO}_aln_edit.fa

#TRIM

#java -jar /data/home/btx494/Programmes/BMGE-1.12/BMGE.jar -i gene_trees/${ORTHO}_aln.fa -t AA -of gene_trees/${ORTHO}_alntrimmed.fa
#java -jar /data/home/btx494/Programmes/BMGE-1.12/BMGE.jar -i gene_trees/${ORTHO}_aln_edit.fa -t AA -o gene_trees/${ORTHO}_aln_edit_trimmed.phy

#MAKE GENE TREES

module load raxml

raxmlHPC-SSE3 	-f a \
		-n ${ORTHO} \
		-s gene_trees/${ORTHO}_aln_edit_trimmed.phy \
		-N 1000 \
		-m PROTGAMMAAUTO \
		-p 12345 \
		-x 12345 \
		-w /data/SBCS-BuggsLab/RowenaHill/fus_comparison/phylogenomics/gene_trees/RAxML

#rm gene_trees/${ORTHO}_aln_edit_trimmed.phy.reduced

