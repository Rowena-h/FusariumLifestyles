#!/bin/sh

cd modeltest-ng

NUM=$(ls -d partition.* | wc -l)

qsub -t 1-${NUM} modeltestng.sh
