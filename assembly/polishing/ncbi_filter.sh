#!/bin/sh

STRAINS=$(cat ../strains)

for STRAIN in $STRAINS
do
	
	awk '{print $1}' mito_${STRAIN} > tmp

	#Remove sequences flagged as mitochondrial contaminants by NCBI
	awk 'BEGIN{while((getline<"tmp")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' fusotu${STRAIN}_abyss_pilon_filtered.fa > tmp1 && mv tmp1 fusotu${STRAIN}_abyss_pilon_filtered.fa
	rm tmp

	if [ -f duplicates_${STRAIN} ]
	then

		awk '{print $1}' duplicates_${STRAIN} > tmp		

		#Remove sequences flagged as duplicates
		awk 'BEGIN{while((getline<"tmp")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' fusotu${STRAIN}_abyss_pilon_filtered.fa > tmp1 && mv tmp1 fusotu${STRAIN}_abyss_pilon_filtered.fa
		rm tmp
	
	fi

done
