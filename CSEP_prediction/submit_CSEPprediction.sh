#!/bin/sh
#Script to submit all programmes in the CSEP prediction pipeline

#Directory containing protein fasta files (MUST END IN FORWARD SLASH)
DIR=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/proteins/
#File listing the files to be run
ls ../proteins/*.faa | sed 's#\.\./proteins/##' > todo
SAMPLES=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/CSEP_prediction/todo
#Number of samples
NUM=$(cat $SAMPLES | wc -l)

#Submit all jobs
qsub phobius/phobius.sh $DIR $SAMPLES
qsub prosite/ps_scan.sh $DIR $SAMPLES
qsub predgpi/predgpi.sh $DIR $SAMPLES

for i in $(cat $SAMPLES)
do
	mkdir targetp/tmp_${i}
done

qsub -t 1-${NUM} targetp/targetp.sh $DIR $SAMPLES
qsub tmhmm/tmhmm.sh $DIR $SAMPLES
qsub -t 1-${NUM} effectorp/effectorp.sh $DIR $SAMPLES
qsub nucpred/nucpred.sh $DIR $SAMPLES
qsub signalp/signalp.sh $DIR $SAMPLES
