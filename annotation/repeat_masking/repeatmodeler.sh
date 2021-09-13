#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 5		# Request 5 cores
#$ -l h_rt=24:00:0 	# Request 24 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -j y
#$ -m bea
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../strains)

mkdir fusotu${STRAIN}

cd fusotu${STRAIN}

module load singularity

singularity exec /data/containers/repeatmodeler/repeatmodeler-2.0.1.simg BuildDatabase -name fusotu${STRAIN}_abyss ../../../assembly/polishing/fusotu${STRAIN}_abyss_pilon.fasta

singularity exec /data/containers/repeatmodeler/repeatmodeler-2.0.1.simg RepeatModeler -database fusotu${STRAIN}_abyss -engine ncbi -pa 1 -LTRStruct >& fusotu${STRAIN}_abyss.out 
