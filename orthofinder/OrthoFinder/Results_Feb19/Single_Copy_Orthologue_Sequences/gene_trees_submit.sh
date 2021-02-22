#!/bin/sh

ls *.fa > aln_list

mkdir trees

qsub gene_trees.sh
