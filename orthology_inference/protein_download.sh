#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=1:00:00

#Read ftp link for proteins into variable
LINK=$(sed -n ${SGE_TASK_ID}p fus_ncbi_proteins)

#Download from the ncbi ftp server
wget $LINK

#Read taxon name into variable
FILE=$(echo $LINK | sed 's/^.*\(GCA.*\.faa\).*$/\1/')

#Extract file
gunzip ${FILE}.gz

#Add filename to fasta headers
sed -i "s/>/>${FILE}_/g" $FILE
