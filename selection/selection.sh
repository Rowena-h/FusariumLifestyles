#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1 		# Request 1 cores
#$ -l h_rt=1:0:0 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -t 1-23
#$ -o /dev/null
#$ -j y

ORTHO=$(cat orthogroups_selection_of.csv | sed -n ${SGE_TASK_ID}p)

module load mafft

mafft /data/SBCS-BuggsLab/RowenaHill/fus_comparison/orthofinder/all_fus/OrthoFinder/Results_Oct21/Orthogroup_Sequences/${ORTHO}.fa  > alignments/${ORTHO}_aln.fa

#Remove jgi prefixes
sed -i 's/jgi|*|.*|//' alignments/${ORTHO}_aln.fa

#Make list of sequences
grep ">" alignments/${ORTHO}_aln.fa | sed 's/>//g' > alignments/${ORTHO}_seqlist

module load anaconda3
conda activate biopython

#Get nucleotides for genbank proteins 
python pull_nucleotides.py alignments/${ORTHO}_seqlist

#Get nucleotides for own proteins
awk -F '>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' alignments/${ORTHO}_seqlist own_fus_transcripts.fa >> alignments/${ORTHO}_seqlist_nucl.fa

perl /data/home/btx494/Programmes/pal2nal.v14/pal2nal.pl alignments/${ORTHO}_aln.fa alignments/${ORTHO}_seqlist_nucl.fa -output paml -nogap > alignments/for_paml/${ORTHO}.codon
perl /data/home/btx494/Programmes/pal2nal.v14/pal2nal.pl alignments/${ORTHO}_aln.fa alignments/${ORTHO}_seqlist_nucl.fa -output fasta -nogap > alignments/${ORTHO}_aln_nuc.fa

#Remove protein name
sed -i 's/\.faa.*//g' alignments/${ORTHO}_aln_nuc.fa

#Check files
if [ `head -1 alignments/for_paml/${ORTHO}.codon | awk '{print int($1)}'` != `grep '>' alignments/${ORTHO}_aln.fa | wc -l` ]
then
        echo ${ORTHO} >> errors.txt
fi
