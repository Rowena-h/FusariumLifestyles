#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 core
#$ -l h_rt=01:00:00     # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

ORTHOS=$(tail -n +2 ../sortadate/dating_orthogroups | sed 's/^.*OG000/OG000/' | sed 's/\.raxml\.bestTree_rooted.*//' | sed 's#^#\.\./\.\./phylogenomics/gene_trees/#' | sed 's/$/_aln_edit_trimmed\.phy/')

module load anaconda3
conda activate AMAS

#Concatenate top ten clock-like genes according to SortaDate
AMAS.py concat -f phylip -d aa -i ${ORTHOS} -p fus_proteins_dating10_partition.txt -t fus_proteins_dating10_mcmctree.phy -u phylip

sed -i 's/[[:space:]]/  /' fus_proteins_dating10_mcmctree.phy
sed -e 's/^\(.\{13\}\).*\( .*\)$/\1 \2/' fus_proteins_dating10_mcmctree.phy > fus_proteins_dating10_mcmctree_short.phy

module load R

#Format rooted species tree for MCMCTree
Rscript blank_topology.r ../fus_proteins_62T_iqtree_genepart.contree_rooted ./

sed -i '1s/^/62 1\n/' fus_proteins_62T_iqtree_genepart.contree_rooted_blank
sed -i 's/Root//' fus_proteins_62T_iqtree_genepart.contree_rooted_blank
sed "s/;/\'>0.9<1.35\';/" fus_proteins_62T_iqtree_genepart.contree_rooted_blank > fus_proteins_62T_iqtree_genepart.contree_dating
sed -i "s/GCA_013266205)))/GCA_013266205)))'>0.5<0.9'/" fus_proteins_62T_iqtree_genepart.treefile_dating

module load anaconda3
conda activate paml

#Run MCMCTree
mcmctree mcmctree_step1.ctl

sed -i 's/aaRatefile =/aaRatefile = wag.dat/' tmp*.ctl
sed -i 's/model = 0/model = 2/' tmp*.ctl

for i in tmp*.ctl
do
	echo -e "fix_alpha = 0\nalpha = 0.5\nncatG = 4" >> ${i}
	codeml ${i}
done	

mv rst2 in.BV
