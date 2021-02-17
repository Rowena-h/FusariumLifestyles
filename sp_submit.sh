#!/bin/sh

DIR=$1
SAMPLES=$2

qsub phobius/phobius.sh $DIR $SAMPLES
qsub prosite/ps_scan.sh $DIR $SAMPLES
qsub targetp/targetp.sh $DIR $SAMPLES
qsub tmhmm/tmhmm.sh $DIR $SAMPLES
qsub /data/home/btx494/Programmes/nucpred-1.1/nucpred.sh $DIR $SAMPLES
qsub /data/home/btx494/Programmes/signalp-5.0b/bin/signalp.sh $DIR $SAMPLES
