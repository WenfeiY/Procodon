#!/bin/bash

prefix=$1
cov=$2

root_dir=$prefix
asm_root_dir=$root_dir/spades.assemble/$cov
mkdir -p $prefix.contigs/$cov
for line in $(cat all_gene_list.txt); do
	if [ -e "$asm_root_dir/$line/contigs.fixed.fasta" ]; then
		cp $asm_root_dir/$line/contigs.fixed.fasta $prefix.contigs/$cov/$line.contigs.fasta &
	elif [ -e "$asm_root_dir/$line/contigs.fasta" ]; then
		cp $asm_root_dir/$line/contigs.fasta $prefix.contigs/$cov/$line.contigs.fasta &
	fi
done
