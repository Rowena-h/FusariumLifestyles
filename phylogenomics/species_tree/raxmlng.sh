#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 20            # Request 1 cores
#$ -l h_rt=240:00:00      # Request 1 hour runtime
#$ -l h_vmem=2G         # Request 1GB RAM
#$ -j y
#$ -m bea

#Combine ModelTest-NG partition results
NUM=$(ls -dir modeltest-ng/partition.* | wc -l)

rm fus_proteins_62T_raxmlpartition.txt

for i in $(seq 1 $NUM)
do
	cat modeltest-ng/partition.${i}/fus_proteins_62T_concat.phy.part.aic >> fus_proteins_62T_raxmlpartition.txt
done

module load anaconda3
conda activate raxml-ng

raxml-ng --parse \
	 --msa fus_proteins_62T_concat.phy \
	 --model fus_proteins_62T_raxmlpartition.txt \
	 --prefix raxml-ng/fus_proteins_62T

raxml-ng --all \
	 --msa fus_proteins_62T_concat.phy \
	 --model fus_proteins_62T_raxmlpartition.txt \
	 --prefix raxml-ng/fus_proteins_62T \
	 --seed 2 \
	 --threads ${NSLOTS} \
	 --bs-trees 100

raxml-ng --bsconverge \
	 --bs-trees raxml-ng/fus_proteins_62T.raxml.bootstraps \
	 --prefix raxml-ng/fus_proteins_62T_convergence_test \
	 --seed 2
