#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=1:00:00               #genomes to download

LINK=$(cat links)

for i in $LINK
do
	wget $i
done 

