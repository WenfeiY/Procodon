#!/bin/bash
prefix=$1
cov=30x

#  Align contigs with protein sequences
bash Attached-11-TBLASTN_cov.sh $prefix $cov
