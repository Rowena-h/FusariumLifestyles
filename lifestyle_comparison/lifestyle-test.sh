#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1           # Request 48 cores
#$ -l h_rt=00:30:00     # Request 72 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y

module load anaconda3
conda activate lifestyle-test

module load R

mkdir CSEPs orthogroups

python run_edited.py 	-i lifestyle-test-orthogroups.csv \
			-t species_tree_ingroup.tre \
			--colors endophyte:#009E73,insectmutualist:#56B4E9,plantpathogen:#696969,saprotroph:#0072B2,plantassociate:#9AE324,mycoparasite:#D55E00 \
			-o orthogroups

python run_edited.py  	-i lifestyle-test-CSEPs.csv \
                	-t species_tree_ingroup.tre \
                	--colors endophyte:#009E73,insectmutualist:#56B4E9,plantpathogen:#696969,saprotroph:#0072B2,plantassociate:#9AE324,mycoparasite:#D55E00 \
	                -o CSEPs
