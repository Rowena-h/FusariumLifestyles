#!/bin/sh

cd species_tree

mkdir astral

qsub astral.sh
qsub astral-pro.sh
