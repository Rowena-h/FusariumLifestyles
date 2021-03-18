#!/bin/sh

ls -1 ../orthology_inference/OrthoFinder/Results_Mar01/Single_Copy_Orthologue_Sequences/*.fa | sed s/^.*\\/\// > aln_list

NUM=$(cat aln_list | wc -l)

mkdir gene_trees
mkdir gene_trees/RAxML

qsub -t 1-${NUM} RAxML_genetrees.sh
