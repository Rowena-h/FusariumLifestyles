#!/bin/sh

cat gbff_files/GCA* > gbff_files/concat.gbff

rm errors.txt
rm own_fus_transcripts.fa

for i in 1 3 5 6 7
do
	cat ../../genome_assemblies/Fusarium_OTU${i}/maker/round3/fusotu${i}_abyss50_rnd3.maker.output/fusotu${i}_abyss50_rnd3.all.maker.transcripts.fasta | sed "s/>/>fusotu${i}_abyss50_rnd3.all.maker.proteins.faa_/g" >> own_fus_transcripts.fa
done

sed -i 's/ .*//' own_fus_transcripts.fa

cat own_fus_transcripts.fa Ilysp1_GeneCatalog_transcripts_20121116.nt.fasta | sed 's/.*|/\>Ilysp1_GeneCatalog_proteins_20121116.faa_/' > tmp && mv tmp own_fus_transcripts.fa

mkdir alignments
mkdir alignments/for_paml

dos2unix orthogroups_selection_of.csv

qsub selection.sh
