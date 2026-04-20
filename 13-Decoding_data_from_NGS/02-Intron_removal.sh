#!/bin/bash
prefix=$1
cov=30x

#  Remove intron from assembled contigs
bash Attached-08-remove_intron_cov.sh $prefix $cov
