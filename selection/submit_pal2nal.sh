#!/bin/sh

#Number of samples
NUM_ORTHO=$(ls -1 ../phylogenomics/gene_trees/*_aln.fa | wc -l)

cat gbff_files/GCA* > gbff_files/concat.gbff

rm gbff_files/own_fus_transcripts.fa
rm pal2nal_check

for i in 1 3 5 6 7
do
	cat ../annotation/maker/round3/fusotu${i}_gag2/fusotu${i}.mrna.fasta | sed 's/ .*//' | sed "s/>/>fusotu${i}.proteins.faa_/g" >> gbff_files/own_fus_transcripts.fa
done

sed -i 's/ .*//' gbff_files/own_fus_transcripts.fa

cat gbff_files/own_fus_transcripts.fa gbff_files/Ilysp1_GeneCatalog_transcripts_20121116.nt.fasta | sed 's/.*|/\>Ilysp1_GeneCatalog_proteins_20121116.faa_/' > tmp && mv tmp gbff_files/own_fus_transcripts.fa

mkdir alignments alignments/codon trees

qsub -t 1-${NUM_ORTHO} pal2nal.sh
