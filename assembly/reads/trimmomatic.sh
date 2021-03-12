#!/bin/sh
#$ -cwd           
#$ -pe smp 4      
#$ -l h_rt=01:00:00 
#$ -l h_vmem=2G   
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load trimmomatic

java -jar /share/apps/centos7/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads ${NSLOTS} -trimlog fusotu${STRAIN}-trimmomatic.log FUS_OTU${STRAIN}_1.fastq.gz FUS_OTU${STRAIN}_2.fastq.gz FUS_OTU${STRAIN}_1_trimmedpaired.fastq.gz FUS_OTU${STRAIN}_1_trimmedunpaired.fastq.gz FUS_OTU${STRAIN}_2_trimmedpaired.fastq.gz ${STRAIN}_2_trimmedunpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
