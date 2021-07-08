#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 cores
#$ -l h_rt=01:00:00     # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

ORTHOS=$(cat ../sortadate/dating_orthogroups | sed 's/^.*OG000/OG000/' | sed 's/_rooted.*//' | sed 's#^#\.\./alignments/codon/#' | sed 's/$/_aln_nuc\.fa/') 

mkdir clock_testing clock_testing/independent clock_testing/correlated

module load anaconda3
conda activate AMAS

AMAS.py concat -f phylip -d dna -i ${ORTHOS} -p fus_nucl_dating10_partition.txt -t fus_nucl_dating10_mcmctree.phy -u phylip

sed -i 's/[[:space:]]/  /' fus_nucl_dating10_mcmctree.phy
sed -e 's/^\(.\{13\}\).*\( .*\)$/\1 \2/' fus_nucl_dating10_mcmctree.phy > fus_nucl_dating10_mcmctree_short.phy

module load R

Rscript blank_topology.r ../fus_proteins_62T_iqtree_genepart.contree_rooted ./

sed -i '1s/^/62 1\n/' fus_proteins_62T_iqtree_genepart.contree_rooted_blank
sed -i 's/Root//' fus_proteins_62T_iqtree_genepart.contree_rooted_blank
sed "s/;/\'B\(0.999,1.001\)\';/" fus_proteins_62T_iqtree_genepart.contree_rooted_blank > clock_testing/fus_proteins_62T_iqtree_genepart.contree_clocktest

cp mcmctree_clocktest.ctl clock_testing/independent/mcmctree_clocktest_independent.ctl
cp mcmctree_clocktest.ctl clock_testing/correlated/mcmctree_clocktest_correlated.ctl

qsub mcmctree_clocktest_independent.sh mcmctree_clocktest_correlated.sh
