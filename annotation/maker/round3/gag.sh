#!/bin/sh
#$ -cwd                 # Set the working directory for the job to the current directory
#$ -pe smp 1            # Request 1 core
#$ -l h_rt=1:00:0       # Request 1 hour runtime
#$ -l h_vmem=1G         # Request 1GB RAM
#$ -j y
#$ -t 3

STRAIN=$(sed -n ${SGE_TASK_ID}p ../../strains)

python2.7 ~/Programmes/genomeannotation-GAG-997e384/gag.py \
-f ../../../assembly/polishing/fusotu${STRAIN}_abyss_pilon_filtered.fa \
-g fusotu${STRAIN}_abyss_rnd3.maker.output/fusotu${STRAIN}_abyss_rnd3.all.maker.gff \
-ris 10 \
--fix_terminal_ns \
--fix_start_stop \
-o fusotu${STRAIN}_gag

cd fusotu${STRAIN}_gag

rename genome fusotu${STRAIN} *
mv fusotu${STRAIN}.fasta fusotu${STRAIN}.fsa
sed -i 's/protein|//' fusotu${STRAIN}.proteins.fasta
