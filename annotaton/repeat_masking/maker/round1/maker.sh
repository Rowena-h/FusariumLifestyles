#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1	# Request 48 cores
#$ -l h_rt=120:00:0 	# Request 72 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -t 1-5

module load maker

maker -base fusotu${STRAIN}_abyss_rnd1 -RM_off fusotu${STRAIN}_abyss_round1_maker_opts.ctl
	
gff3_merge -s -d fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1_master_datastore_index.log > fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1.all.maker.gff
fasta_merge -d fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1_master_datastore_index.log
gff3_merge -n -s -d fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1_master_datastore_index.log > fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1.all.maker.noseq.gff
