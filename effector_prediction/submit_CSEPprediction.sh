#!/bin/sh

#Directory containing protein fasta files (MUST END IN FORWARD SLASH)
DIR=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/orthology_inference/
#File listing the protein file names
SAMPLES=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/effector_prediction/todo
#Number of samples
NUM=$(cat $SAMPLES | wc -l)

#qsub phobius/phobius.sh $DIR $SAMPLES
#qsub prosite/ps_scan.sh $DIR $SAMPLES

#for i in $(cat $SAMPLES)
#do
#	mkdir targetp/tmp_${i}
#done

#qsub -t 1-${NUM} targetp/targetp.sh $DIR $SAMPLES
#qsub tmhmm/tmhmm.sh $DIR $SAMPLES
qsub -t 1-${NUM} effectorp/effectorp.sh $DIR $SAMPLES
#qsub nucpred/nucpred.sh $DIR $SAMPLES
#qsub signalp/signalp.sh $DIR $SAMPLES
