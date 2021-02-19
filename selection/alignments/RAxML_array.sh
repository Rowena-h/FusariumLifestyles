#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1		# Request 1 cores
#$ -l h_rt=1:00:00 	# Request 1 hour runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -t 18
#$ -j y

ORTHO=$(cat ../orthogroups_selection.csv | sed -n ${SGE_TASK_ID}p)

perl /data/home/btx494/Scripts/Fasta2Phylip.pl ${ORTHO}_aln_nuc.fa ${ORTHO}_aln_nuc.phy

LENGTH=$(head -1 ${ORTHO}_aln_nuc.phy | awk '{print int($2)}')

echo "DNA, codon1 = 1-`expr $LENGTH - 2`\3" > ${ORTHO}.partition
echo "DNA, codon2 = 2-`expr $LENGTH - 1 `\3" >> ${ORTHO}.partition
echo "DNA, codon3 = 3-${LENGTH}\3" >> ${ORTHO}.partition

module load raxml

raxmlHPC-SSE3   -f a \
                -n ${ORTHO} \
                -s ${ORTHO}_aln_nuc.phy \
		 -q ${ORTHO}.partition \
                -N 1000 \
                -m GTRGAMMA \
                -p 12345 \
                -x 12345

module load R

Rscript ~/Scripts/reroot.r RAxML_bipartitions.${ORTHO} "Ilysp1_GeneCatalog_proteins_20121116"
