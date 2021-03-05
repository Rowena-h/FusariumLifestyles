#!/bin/sh

#Directory containing protein fasta files (MUST END IN FORWARD SLASH)
DIR=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/orthofinder/blah/
#File listing the protein file names
SAMPLES=/data/SBCS-BuggsLab/RowenaHill/fus_comparison/todo2

qsub phobius/phobius.sh $DIR $SAMPLES
qsub prosite/ps_scan.sh $DIR $SAMPLES
qsub targetp/targetp.sh $DIR $SAMPLES
qsub tmhmm/tmhmm.sh $DIR $SAMPLES
qsub /data/home/btx494/Programmes/nucpred-1.1/nucpred.sh $DIR $SAMPLES
qsub /data/home/btx494/Programmes/signalp-5.0b/bin/signalp.sh $DIR $SAMPLES
