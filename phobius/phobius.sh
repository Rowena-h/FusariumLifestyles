#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1    	# Request 48 cores
#$ -l h_rt=24:00:0     # Request 72 hours runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -m bea

DIR=$1
SAMPLES=$(cat $2)
RUN_DIR=/data/SBCS-BuggsLab/RowenaHill/fus_comparison

for i in $SAMPLES
do
	/data/home/btx494/Programmes/phobius/phobius.pl ${DIR}${i} > ${RUN_DIR}/phobius/${i}_phobius
	#Filter for list of SPs	
	grep -i -B 1 "SIGNAL" ${RUN_DIR}/phobius/${i}_phobius | grep "ID" | awk '{print $2}' > ${RUN_DIR}/phobius/${i}_phobius_SPlist
	#Filter for list of >1 TMs
	#Combine Phobius output into one line per gene
	awk 'NR==1{printf $0" ";next}{printf /^ID/ ? "\n"$0" " : $0}' ${RUN_DIR}/phobius/${i}_phobius > tmp
	#Only print rows with more than one TM
	cat tmp | grep -o -n "TRANSMEM" | cut -d : -f 1 | uniq -c | awk '$1>1 {print $2}' | sed '/[^0-9]/d;s/.$/&p/' | sed -nf - tmp | awk '{print $2}' > ${RUN_DIR}/phobius/${i}_phobius_TMlist
	rm tmp
done
