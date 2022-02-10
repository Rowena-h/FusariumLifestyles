#!/bin/sh
#Script to submit blast of CSEPs against the PHI-base database

module load blast+

makeblastdb -dbtype prot -in blastp/phi-base_current.fas

#Number of samples
NUM=$(ls ../proteins/*.faa | sed 's#\.\./proteins/##' | wc -l)

cd blastp

qsub -t 1-${NUM} blastp.sh

