#!/bin/sh

ls -1 ../orthology_inference/OrthoFinder/Results_Mar01/Single_Copy_Orthologue_Sequences/*.fa | sed s/^.*\\/\// > aln_list

mkdir gene_trees
mkdir gene_trees/RAxML

qsub RAxML_genetrees.sh
