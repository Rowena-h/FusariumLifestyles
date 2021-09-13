#!/bin/sh
#Script to submit assembly assessment jobs

qsub quast.sh
qsub busco.sh
qsub blast.sh