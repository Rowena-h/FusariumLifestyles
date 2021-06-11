#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1           # Request 48 cores
#$ -l h_rt=00:30:00     # Request 72 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

module load anaconda3
conda activate lifestyle-test

module load R

mkdir effectors

python run.py 	-i lifestyle-test-effectors.csv \
		-t species_tree_ingroup.tre \
		--colors endophyte:#009E73,insectassociated:#56B4E9,plantpathogen:#696969,saprotroph:#0072B2,plantassociated:#000000,mycoparasite:#D55E00 \
		-o effectors

#python run.py -i raxml\fus_orthogroups\lifestyle-test-orthogroups.csv -t 
#raxml\fus_orthogroups\RAxML_bipartitions.fus_proteins_27T_bsadd_rooted_ingroup --colors 
#endophyte:#009E73,insectmutualist:#56B4E9,plantpathogen:#696969,saprotroph:#0072B2 -o "D:\Documents\Bio 
#Programmes\Effect-Of-Biological-Categories-On-Genomes-Composition-master\Effect-Of-Biological-Categories-On-Genomes-Composition-master\raxml\fus_orthogroups"
