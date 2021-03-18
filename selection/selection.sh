#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1 		# Request 1 cores
#$ -l h_rt=1:0:0 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y

ORTHO=$(ls -1 ../phylogenomics/gene_trees/*_aln.fa | sed 's#\.\./phylogenomics/gene_trees/##' | sed 's/_aln\.fa//' | sed -n ${SGE_TASK_ID}p)

#Remove jgi prefixes
sed 's/jgi|*|.*|//' ../phylogenomics/gene_trees/${ORTHO}_alntrimmed.fa > alignments/${ORTHO}_aln.fa

#Make list of sequences
grep ">" alignments/${ORTHO}_aln.fa | sed 's/>//g' > alignments/${ORTHO}_seqlist

module load anaconda3
conda activate biopython

#Get nucleotides for genbank proteins 
python pull_nucleotides.py alignments/${ORTHO}_seqlist

#Get nucleotides for own proteins
awk -F '>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' alignments/${ORTHO}_seqlist own_fus_transcripts.fa >> alignments/${ORTHO}_seqlist_nucl.fa

perl /data/home/btx494/Programmes/pal2nal.v14/pal2nal.pl alignments/${ORTHO}_aln.fa alignments/${ORTHO}_seqlist_nucl.fa -output fasta -nogap > alignments/${ORTHO}_aln_nuc.fa

#Remove protein name
sed -i 's/\.faa.*//g' alignments/${ORTHO}_aln_nuc.fa

module load R

#Reroot gene trees
Rscript reroot.r RAxML_bipartitions.${ORTHO} "Ilysp1_GeneCatalog_proteins_20121116" /data/SBCS-BuggsLab/RowenaHill/fus_comparison/selection/alignments/
