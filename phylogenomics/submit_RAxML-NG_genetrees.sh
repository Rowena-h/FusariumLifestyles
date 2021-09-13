#!/bin/sh

NUM=$(cat aln_list | wc -l)

mkdir gene_trees/RAxML-NG

qsub -t 1-${NUM} RAxMLNG_genetrees.sh
