#!/bin/sh

#Number of samples
NUM=$(ls ../proteins/*.faa | sed 's#\.\./proteins/##' | wc -l)

qsub -t 1-${NUM} run_dbcan.sh/run_dbcan.sh
