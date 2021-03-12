#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1	# Request 48 cores
#$ -l h_rt=120:00:0 	# Request 72 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../../../strains)

module load maker
export AUGUSTUS_CONFIG_PATH=/data/SBCS-BuggsLab/RowenaHill/genome_assemblies/augustus_config/config

# transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' ../round1/fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1.all.maker.noseq.gff > ../round1/fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1.all.maker.est2genome.gff
# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' ../round1/fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1.all.maker.noseq.gff > ../round1/fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1.all.maker.protein2genome.gff

maker -base fusotu${STRAIN}_abyss_rnd2 -RM_off fusotu${STRAIN}_abyss_round2_maker_opts.ctl
   
gff3_merge -s -d fusotu${STRAIN}_abyss_rnd2.maker.output/fusotu${STRAIN}_abyss_rnd2_master_datastore_index.log > fusotu${STRAIN}_abyss_rnd2.maker.output/fusotu${STRAIN}_abyss_rnd2.all.maker.gff
fasta_merge -d fusotu${STRAIN}_abyss_rnd2.maker.output/fusotu${STRAIN}_abyss_rnd2_master_datastore_index.log
gff3_merge -n -s -d fusotu${STRAIN}_abyss_rnd2.maker.output/fusotu${STRAIN}_abyss_rnd2_master_datastore_index.log > fusotu${STRAIN}_abyss_rnd2.maker.output/fusotu${STRAIN}_abyss_rnd2.all.maker.noseq.gff
