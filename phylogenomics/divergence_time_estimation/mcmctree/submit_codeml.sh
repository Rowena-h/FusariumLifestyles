#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 cores
#$ -l h_rt=00:30:00     # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

#ORTHOS=$(cat ../SortaDate/dating_orthogroups | sed 's/^.*OG000/OG000/' | sed 's/_rooted.*//' | sed 's#^#\.\./\.\./gene_trees/#' | sed 's/$/_aln_edit_trimmed\.phy/')

#rm fus_proteins_dating10_mcmctree.phy

#for i in $ORTHOS
#do
#    cat ${i} >> fus_proteins_dating10_mcmctree.phy
#done

#sed -i 's/[[:space:]]/  /' fus_proteins_dating10_mcmctree.phy
#sed -e 's/^\(.\{13\}\).*\( .*\)$/\1 \2/' fus_proteins_dating10_mcmctree.phy > fus_proteins_dating10_mcmctree_short.phy

module load R

#Rscript ../../../selection/reroot.r ../../species_tree/iqtree/gene_partitions/fus_proteins_62T_iqtree_genepart.treefile Ilysp1_GeneCatalog_proteins_20121116 ./
#Rscript blank_topology.r fus_proteins_62T_iqtree_genepart.treefile_rooted ./

#sed -i '1s/^/62 1\n/' fus_proteins_62T_iqtree_genepart.treefile_rooted_blank
#sed -i 's/Root//' fus_proteins_62T_iqtree_genepart.treefile_rooted_blank

module load anaconda3
conda activate paml

mcmctree mcmctree_step1.ctl

wait

mkdir codeml

cp tmp00* codeml

for i in tmp00*.trees
do
	Rscript ../../../selection/reroot.r codeml/${i} "Ilysp1_GeneCa" codeml/
done	

rename trees_rooted trees codeml/tmp00*


for i in tmp00*.ctl
do
	echo "clock = 1" >> codeml/${i}
	sed -i 's/getSE = 2/getSE = 0/' codeml/${i}
	codeml codeml/${i}
done

sed -i 's/aaRatefile =/aaRatefile = wag.dat/' tmp00*.ctl
sed -i 's/model = 0/model = 2/' tmp00*.ctl

for i in tmp00*.ctl
do
	echo -e "fix_alpha = 0\nalpha = 0.5\nncatG = 4" >> ${i}
	codeml $i

done	

wait

mv rst2 in.BV

#module load anaconda3
#source activate AMAS

#AMAS.py concat -f phylip -d aa -i ${ORTHOS} -p fus_proteins_dating10_partition.txt -t fus_proteins_dating10.phy -u phylip
