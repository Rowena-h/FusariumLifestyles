#!/bin/sh

module load bedtools

STRAINS=$(cat ../strains)

for STRAIN in $STRAINS
do

	STRAIN=5	
	
	awk '{print $1}' mito_remove_${STRAIN} > tmp

	#Remove sequences flagged as mitochondrial contaminants by NCBI
	awk 'BEGIN{while((getline<"tmp")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' fusotu${STRAIN}_abyss_pilon_filtered.fa > tmp1 && mv tmp1 fusotu${STRAIN}_abyss_pilon_filtered.fa
	rm tmp

	#Trim regions flagged as mitochondrial contaminants by NCBI
	bedtools maskfasta -fi fusotu${STRAIN}_abyss_pilon_filtered.fa -bed mito_trim_${STRAIN}.bed -fo tmp.fa -mc X
	sed 's/X//g' tmp.fa > fusotu${STRAIN}_abyss_pilon_filtered.fa
	rm tmp.fa

	if [ -f duplicates_${STRAIN} ]
	then

		awk '{print $1}' duplicates_${STRAIN} > tmp		

		#Remove sequences flagged as duplicates
		awk 'BEGIN{while((getline<"tmp")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' fusotu${STRAIN}_abyss_pilon_filtered.fa > tmp1 && mv tmp1 fusotu${STRAIN}_abyss_pilon_filtered.fa
		rm tmp
	
	fi

done
