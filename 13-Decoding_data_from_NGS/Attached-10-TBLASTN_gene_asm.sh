#!/bin/bash

#  Read prefix, coverage and gene name
prefix=$1
gene=$2

#  The directory of protein sequences from yeast genome
protein_root_dir="./S.cerevisiae.Protein"

#  Default directories of results
root_dir=$prefix
TBLASTN_asm_dir=$root_dir/TBLASTN_asm_out

#  File path of contigs
contig=$prefix.contigs/$gene.contigs.fasta

#  FASTA file path of protein sequence
protein_seq_path=$protein_root_dir/$gene.Protein.fa

#  TBLASTN results of protein sequences aligned with contigs
TBLASTN_asm_out=$TBLASTN_asm_dir/$gene.tblastn.out

#  TBLASTN (protein sequences to pair-ended reads)
tblastn -query $protein_seq_path \
    -subject $contig \
    -outfmt 5 \
    -evalue 0.1 \
    -matrix PAM30 \
    -soft_masking true \
    -gapopen 15 \
    -gapextend 3 \
    > $TBLASTN_asm_out
