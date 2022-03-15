#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 1 core
#$ -l h_rt=1:00:00 	# Request 30 minutes runtime
#$ -l h_vmem=3G   	# Request 1GB RAM
#$ -j y

mkdir OF_gene_trees

cp ../../orthology_inference/OrthoFinder/Results_Jan28/Resolved_Gene_Trees/*_tree.txt OF_gene_trees

find OF_gene_trees/ -type f -exec sed -i 's/proteins_[^:]*:/proteins:/g' {} \;
find OF_gene_trees/ -type f -exec sed -i 's/protein_[^:]*:/protein:/g' {} \;
find OF_gene_trees/ -type f -exec sed -i -e '$a\' {} \;

TREES=$(ls ../../orthology_inference/OrthoFinder/Results_Jan28/Resolved_Gene_Trees/*_tree.txt | wc -l)

cat OF_gene_trees/*_tree.txt > astral/${TREES}_fusortho_multicopy_OF_trees.tre

~/Programmes/ASTER-master/bin/astral-pro 	-t ${NSLOTS} \
						-o astral/fus_proteins_62T_astralpro_multicopy.tre astral/${TREES}_fusortho_multicopy_OF_trees.tre
