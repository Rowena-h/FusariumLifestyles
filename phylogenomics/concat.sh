#!/bin/sh

mkdir modeltest-ng

module load anaconda3
source activate AMAS

#Concatenate gene alignments
AMAS.py concat -f phylip -d aa -i gene_trees/*_aln_edit_trimmed.phy -p fus_proteins_62T_partition.txt -t fus_proteins_62T_concat.phy -u phylip

#Fix format
sed -i 's/^/PROT, /' fus_proteins_62T_partition.txt

#Split partitions for ModelTest-NG
split -l 100 --numeric-suffixes=1 fus_proteins_62T_partition.txt fus_proteins_62T_partition.num
rename num0 num *.num*

#Prepare partitioned ModelTest-NG folders
NUM=$(ls fus_proteins_62T_partition.num* | wc -l)

for i in $(seq 1 $NUM)
do
	mkdir modeltest-ng/partition.${i}
	mv fus_proteins_62T_partition.num${i} modeltest-ng/partition.${i}
	cp fus_proteins_62T_concat.phy modeltest-ng/partition.${i}
done

#Prepare IQ-TREE partition file
mv fus_proteins_62T_partition.txt species_tree/fus_proteins_62T_iqtreepartition.txt
mv fus_proteins_62T_concat.phy species_tree/

deactivate AMAS
