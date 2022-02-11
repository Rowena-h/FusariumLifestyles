#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 core
#$ -l h_rt=00:30:00 	# Request 30 minutes runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y
#$ -o /dev/null

ORTHO=$(cat aln_list | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)
ORTHOFINDER_DIR=../orthology_inference/OrthoFinder/Results_Oct22/Single_Copy_Orthologue_Sequences

#Align single copy orthogroups

module load mafft

mafft ${ORTHOFINDER_DIR}/${ORTHO}.fa > gene_trees/${ORTHO}_aln.fa

#Fix names

sed 's/.faa.*$//' gene_trees/${ORTHO}_aln.fa > gene_trees/${ORTHO}_aln_edit.fa

#Trim alignments

module load anaconda3
conda activate trimal

trimal -in gene_trees/${ORTHO}_aln.fa -fasta -gappyout > gene_trees/${ORTHO}_alntrimmed.fa
trimal -in gene_trees/${ORTHO}_aln_edit.fa -phylip -gappyout > gene_trees/${ORTHO}_aln_edit_trimmed.phy

java -jar /data/home/btx494/Programmes/BMGE-1.12/BMGE.jar 	-i gene_trees/${ORTHO}_aln.fa \
								-t AA -of gene_trees_bmge/${ORTHO}_alntrimmed.fa
java -jar /data/home/btx494/Programmes/BMGE-1.12/BMGE.jar 	-i gene_trees/${ORTHO}_aln_edit.fa \
								-t AA -o gene_trees_bmge/${ORTHO}_aln_edit_trimmed.phy
