#!/bin/sh
#$ -cwd           	# Set the working directory for the job to the current directory
#$ -pe smp 1	    # Request 1 core
#$ -l h_rt=12:00:0 	# Request 12 hours runtime
#$ -l h_vmem=1G   	# Request 1GB RAM
#$ -m bea
#$ -t 1-5

STRAIN=$(sed -n ${SGE_TASK_ID}p ../../../strains)

module load anaconda3
conda activate snap
module load maker

mkdir ${STRAIN}
cd ${STRAIN}	

# export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -x 0.25 -l 50 -d ../../../round1/fusotu${STRAIN}_abyss_rnd1.maker.output/fusotu${STRAIN}_abyss_rnd1_master_datastore_index.log 
rename genome fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25 *

# gather some stats and validate
fathom fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25.ann fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1 
fathom fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25.ann fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25.dna -validate > validate.log 2>&1

# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25.ann fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorise.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
	
# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
	
# assemble the HMM
hmm-assembler.pl fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25 params > fusotu${STRAIN}_abyss_rnd1.zff.length50_aed0.25.hmm
