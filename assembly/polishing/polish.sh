#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 4		# Request 1 core
#$ -l h_rt=48:0:0 	# Request 24 hour runtime
#$ -l h_vmem=5G   	# Request 1GB RAM
#$ -m bea
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

module load bwa
module load samtools

for ASSEMBLER in abyss megahit spades
do
	#Index assembly	
	bwa index ../denovo_assembly/${ASSEMBLER}/fusotu${STRAIN}/fusotu${STRAIN}-contigs.fa
	#Align reads to assembly and sort by readname
	bwa mem ../denovo_assembly/${ASSEMBLER}/fusotu${STRAIN}/fusotu${STRAIN}-contigs.fa ../FUS_OTU${STRAIN}_1_trimmedpaired.fastq.gz ../FUS_OTU${STRAIN}_2_trimmedpaired.fastq.gz -t ${NSLOTS} | samtools sort -@ ${NSLOTS} -n -o fusotu${STRAIN}_${ASSEMBLER}_mapped_sorted.bam -
	#Coordinate sort and mark duplicates
	samtools fixmate -m -@ ${NSLOTS} fusotu${STRAIN}_${ASSEMBLER}_mapped_sorted.bam - | samtools sort -@ ${NSLOTS} - | samtools markdup -@ ${NSLOTS} - fusotu${STRAIN}_${ASSEMBLER}_mapped_sorted_dups.bam
	#Calculate statistics
	samtools flagstat fusotu${STRAIN}_${ASSEMBLER}_mapped_sorted_dups.bam > fusotu${STRAIN}_${ASSEMBLER}_mapstats
	
	#Coordinate sort for polishing
	samtools sort -@ ${NSLOTS} fusotu${STRAIN}_${ASSEMBLER}_mapped_sorted.bam -o fusotu${STRAIN}_${ASSEMBLER}_mapped_coordinatesorted.bam
	#Index
	samtools index fusotu${STRAIN}_${ASSEMBLER}_mapped_coordinatesorted.bam

	#Polish with pilon
	java -jar /data/home/btx494/pilon/pilon-1.23.jar --genome ../denovo_assembly/${ASSEMBLER}/fusotu${STRAIN}/fusotu${STRAIN}-contigs.fa --frags fusotu${STRAIN}_${ASSEMBLER}_mapped_coordinatesorted.bam --output fusotu${STRAIN}_${ASSEMBLER}_pilon --changes --fix all --threads ${NSLOTS}
done
