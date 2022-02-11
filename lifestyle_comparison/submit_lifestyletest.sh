#!/bin/sh

module load R/3.6.1

Rscript lifestyle_v_phylogeny.r ../CSEP_CAZyme_prediction/orthogroup-matrices-2022-02-10.RData
