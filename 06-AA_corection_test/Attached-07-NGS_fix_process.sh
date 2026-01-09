#!/bin/bash

#  Read prefix, coverage and gene name
prefix=$1
cov=$2
gene=$3

#  The directory of protein sequences from yeast genome
protein_root_dir="./S.cerevisiae.Protein"

#  Default directories of results
root_dir=$prefix
TBLASTN_dir=$root_dir/TBLASTN_out/$cov
FQ_sub_dir=$root_dir/FQ_sub/$cov
FQ_sub_fix_dir=$root_dir/FQ_sub_fix/$cov

#  Create output directory
mkdir -p $FQ_sub_fix_dir 2>/dev/null

#  FASTA file path of protein sequence
protein_seq_path=$protein_root_dir/$gene.Protein.fa

#  TBLASTN results of reads aligned to protein sequences
TBLASTN1_out=$TBLASTN_dir/$gene.tblastn.1.out
TBLASTN2_out=$TBLASTN_dir/$gene.tblastn.2.out

#  File path of FASTQ subset
FQ1_sub=$FQ_sub_dir/$gene.read1.fastq
FQ2_sub=$FQ_sub_dir/$gene.read2.fastq

#  File path of corrected FASTQ subset
FQ1_sub_fix=$FQ_sub_fix_dir/$gene.read1.fastq
FQ2_sub_fix=$FQ_sub_fix_dir/$gene.read2.fastq

#  Run reads correction program
python Attached-06-NGS_fix.py $protein_seq_path $TBLASTN1_out $FQ1_sub $FQ1_sub_fix
python Attached-06-NGS_fix.py $protein_seq_path $TBLASTN2_out $FQ2_sub $FQ2_sub_fix
