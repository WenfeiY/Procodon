#!/bin/bash

prefix=$1

root_dir=$prefix
asm_root_dir=$root_dir/spades.assemble
mkdir -p $prefix.contigs
for line in $(cat all_gene_list.txt); do
	if [ -e "$asm_root_dir/$line/contigs.fixed.fasta" ]; then
		cp $asm_root_dir/$line/contigs.fixed.fasta $prefix.contigs/$line.contigs.fasta &
	elif [ -e "$asm_root_dir/$line/contigs.fasta" ]; then
		cp $asm_root_dir/$line/contigs.fasta $prefix.contigs/$line.contigs.fasta &
	fi
done
