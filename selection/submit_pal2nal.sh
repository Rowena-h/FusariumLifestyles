#!/bin/sh

#Number of samples
#NUM_ORTHO=$(ls -1 ../phylogenomics/gene_trees/*_aln.fa | sed 's#\.\./phylogenomics/gene_trees/##' | sed 's/_aln\.fa//' | wc -l)
NUM_ORTHO=$(cat torepeat | wc -l)

cat gbff_files/GCA* > gbff_files/concat.gbff

#rm errors.txt
rm own_fus_transcripts.fa
rm pal2nal_check

for i in 1 3 5 6 7
do
	cat ../annotation/maker/round3/fusotu${i}_abyss_rnd3.maker.output/fusotu${i}_abyss_rnd3.all.maker.transcripts.fasta | sed "s/>/>fusotu${i}_abyss_rnd3.all.maker.proteins.faa_/g" >> own_fus_transcripts.fa
done

sed -i 's/ .*//' own_fus_transcripts.fa

cat own_fus_transcripts.fa Ilysp1_GeneCatalog_transcripts_20121116.nt.fasta | sed 's/.*|/\>Ilysp1_GeneCatalog_proteins_20121116.faa_/' > tmp && mv tmp own_fus_transcripts.fa

mkdir alignments
mkdir alignments/codon
mkdir trees
#mkdir alignments/for_paml

#dos2unix orthogroups_selection.csv

qsub -t 1-${NUM_ORTHO} pal2nal.sh
