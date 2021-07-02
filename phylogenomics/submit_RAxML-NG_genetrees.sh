#!/bin/sh

NUM=$(cat aln_list2 | wc -l)

#mkdir gene_trees
#mkdir gene_trees/RAxML-NG

qsub -t 1-${NUM} RAxMLNG_genetrees.sh
