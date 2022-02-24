#!/bin/sh

NUM=$(cat aln_list | wc -l)

mkdir gene_trees/RAxML-NG
mkdir gene_trees_bmge/RAxML-NG

qsub -t 1-${NUM} RAxMLNG_genetrees.sh
