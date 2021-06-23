#!/bin/sh

NUM=$(ls -dir partition.* | wc -l)

qsub -t 1-${NUM} modeltestng.sh
