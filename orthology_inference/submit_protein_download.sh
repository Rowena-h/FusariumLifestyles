#!/bin/bash

module load R 

#Get FTP links for Fusarium taxa off NCBI
Rscript ncbi_ftp_links.r

#Read number of taxa into variable
NUM=$(cat fus_ncbi_proteins | wc -l)

#Submit download
qsub -t 1-${NUM} protein_download.sh

#Copy across own proteins
cp ../annotation/maker/round3/fusotu*_gag/fusotu*.proteins.fasta .
rename fasta faa *

for FILE in $(ls fusotu*)
do
	#Add filename to fasta headers
	sed -i "s/>/>${FILE}_/g" $FILE
done
