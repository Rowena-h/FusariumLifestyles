#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 core
#$ -l h_rt=12:00:00 	# Request 12 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

ORTHO=$(cat ../../phylogenomics/aln_list | sed 's/\.fa//' | sed -n ${SGE_TASK_ID}p)

module load anaconda3
conda activate hyphy-2.5.30

for LIFESTYLE in endophyte insectmutualist mycoparasite plantassociate plantpathogen saprotroph
do

	hyphy contrast-fel 	--alignment ../alignments/codon/${ORTHO}_aln_nuc.fa \
                        	--tree ../trees/${ORTHO}.raxml.bestTree_rooted.${LIFESTYLE} \
				--output contrast-fel/${ORTHO}_Contrast-FEL_${LIFESTYLE}.json \
        	                --branch-set lifestyle
done
