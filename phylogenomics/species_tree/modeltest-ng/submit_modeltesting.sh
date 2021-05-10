#!/bin/sh

#module load anaconda3
#conda activate AMAS

#AMAS.py concat -f phylip -d aa -i ../../gene_trees/*_aln_edit_trimmed.phy -p fus_proteins_62T_partition.txt -t fus_proteins_62T_concat.phy -u phylip

#sed -i 's/^/PROT, /' fus_proteins_62T_partition.txt

split -l 100 --numeric-suffixes=1 fus_proteins_62T_partition.txt fus_proteins_62T_partition.num
rename num0 num *.num*

NUM=$(ls fus_proteins_62T_partition.num* | wc -l)

for i in $(seq 1 $NUM)
do
	mkdir partition.${i}
	mv fus_proteins_62T_partition.num${i} partition.${i}
	cp fus_proteins_62T_concat.phy partition.${i}
done

#qsub -t 1-${NUM} modeltestng.sh

#deactivate AMAS
