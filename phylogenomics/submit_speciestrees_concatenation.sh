#!/bin/sh

cd species_tree

qsub raxmlng.sh
qsub iqtree.sh
