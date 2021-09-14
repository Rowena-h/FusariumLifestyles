#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 core
#$ -l h_rt=00:30:00     # Request 30 minutes runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

mkdir ../trees

module load R

#Reroot gene trees for SortaDate
for ORTHO in $(cat ../../phylogenomics/aln_list | sed 's/\.fa//')
do
	Rscript ../reroot.r ../../phylogenomics/gene_trees/RAxML-NG/${ORTHO}.raxml.bestTree "Ilysp1_GeneCatalog_proteins_20121116" ../trees/
done

#Reroot species tree for SortaDate
Rscript ../reroot.r ../../phylogenomics/species_tree/iqtree/fus_proteins_62T_iqtree_genepart.contree "Ilysp1_GeneCatalog_proteins_20121116" ../

module load anaconda3
conda activate phyx

#Run SortaDate
python get_var_length.py ../trees/ --flend _rooted --outf var --outg Ilysp1_GeneCatalog_proteins_20121116

python get_bp_genetrees.py ../trees/ ../fus_proteins_62T_iqtree_genepart.contree_rooted --flend _rooted --outf bp

python combine_results.py var bp --outf comb

python get_good_genes.py comb --max 10 --outf dating_orthogroups
