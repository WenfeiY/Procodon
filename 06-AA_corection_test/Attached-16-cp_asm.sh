#!/bin/bash

for prefix in Train_Ng_0.05 Train_Ng_0.02 Train_Ng_0.01 Train_Ng_0.005
do
	root_dir=$prefix
	for cov in 30x
	do
		asm_root_dir=$root_dir/spades.assemble/$cov
		mkdir -p $prefix.contigs/$cov
		for line in $(cat Train_gene_list.txt); do
			if [ -e "$asm_root_dir/$line/contigs.fixed.fasta" ]; then
				cp $asm_root_dir/$line/contigs.fixed.fasta $prefix.contigs/$cov/$line.contigs.fasta &
			elif [ -e "$asm_root_dir/$line/contigs.fasta" ]; then
				cp $asm_root_dir/$line/contigs.fasta $prefix.contigs/$cov/$line.contigs.fasta &
			fi
		done
	done
done
