#!/bin/bash

#  Read prefix, coverage and gene name
prefix=$1
cov=$2
gene=$3

#  Default directories of results
root_dir=$prefix
FQ_sub_fix_dir=$root_dir/FQ_sub_fix/$cov
asm_fix_root_dir=$root_dir/spades.assemble.fixed/$cov
spades_dir_fix=$asm_fix_root_dir/$gene

#  Corrected NGS data path
FQ1_sub_fix=$FQ_sub_fix_dir/$gene.read1.fastq
FQ2_sub_fix=$FQ_sub_fix_dir/$gene.read2.fastq

#  Assemble corrected reads
spades.py -o $spades_dir_fix/ -1 $FQ1_sub_fix -2 $FQ2_sub_fix --careful -k 21,33,55 -t 8 
