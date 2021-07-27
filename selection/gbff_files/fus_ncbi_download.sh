#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=1:00:00
#$ -t 1-56		#number of taxa to download

#Read ftp link for gbff into variable
LINK=$(sed -n ${SGE_TASK_ID}p ../../fus_ncbi_genomic)

#Download from the ncbi ftp server
wget $LINK

#Read taxon name into variable
FILE=$(echo $LINK | sed 's/^.*\(GCA.*\.gbff\).*$/\1/')

#Extract file
gunzip ${FILE}.gz
