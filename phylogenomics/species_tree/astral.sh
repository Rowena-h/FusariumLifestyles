#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 core
#$ -l h_rt=0:30:00 	# Request 30 minutes runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

mkdir astral

TREES=$(ls ../gene_trees/RAxML-NG/*bestTree | wc -l)

cat ../gene_trees/RAxML-NG/*bestTree > astral/${TREES}_fusortho_raxmlng_trees.tre

java -jar /data/home/btx494/Programmes/Astral/astral.5.7.3.jar 	-i astral/${TREES}_fusortho_raxmlng_trees.tre \
                                                                -o astral/fus_astral_proteins_62T.tre
