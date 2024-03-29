#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1	    	# Request 1 core
#$ -l h_rt=120:00:0 	# Request 120 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -j y
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../../strains)

module load maker
export AUGUSTUS_CONFIG_PATH=/data/SBCS-BuggsLab/RowenaHill/genome_assemblies/augustus_config/config

maker -base fusotu${STRAIN}_abyss_rnd3 -RM_off fusotu${STRAIN}_abyss_round3_maker_opts.ctl
    
gff3_merge -s -d fusotu${STRAIN}_abyss_rnd3.maker.output/fusotu${STRAIN}_abyss_rnd3_master_datastore_index.log > fusotu${STRAIN}_abyss_rnd3.maker.output/fusotu${STRAIN}_abyss_rnd3.all.maker.gff
fasta_merge -d fusotu${STRAIN}_abyss_rnd3.maker.output/fusotu${STRAIN}_abyss_rnd3_master_datastore_index.log
gff3_merge -n -s -d fusotu${STRAIN}_abyss_rnd3.maker.output/fusotu${STRAIN}_abyss_rnd3_master_datastore_index.log > fusotu${STRAIN}_abyss_rnd3.maker.output/fusotu${STRAIN}_abyss_rnd3.all.maker.noseq.gff
