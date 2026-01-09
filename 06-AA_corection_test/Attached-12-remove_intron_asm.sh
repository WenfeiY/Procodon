#!/bin/bash

#  Read prefix, coverage and gene name
prefix=$1
cov=$2
gene=$3

#  The directory of intron sequences from yeast genome
intron_root_dir="./S.cerevisiae.Intron"

#  Default directories of results
root_dir=$prefix
BLASTN_dir=$root_dir/BLASTN_out/$cov
asm_root_dir=$root_dir/spades.assemble/$cov
spades_dir=$asm_root_dir/$gene

#  FASTA file path of intron sequence
intron_seq_path=$intron_root_dir/$gene.intron.fa
#  BLASTN result path of introns aligned with contigs
asm_BLASTN_intron_out=$BLASTN_dir/$gene.asm.intron.blastn.out

#  BLASTN align introns and contigs
blastn -query $intron_seq_path \
    -subject $spades_dir/contigs.fasta \
    -outfmt 6 \
    -word_size 15 \
    -evalue 1e-5 \
    > $asm_BLASTN_intron_out

#  Default path of assembled contigs
asm_path=$spades_dir/contigs.fasta
#  Default path of output contigs (intron removed)
asm_fix_path=$spades_dir/contigs.fixed.fasta

#  Remove intron from contigs
python Attached-11-remove_intron_from_contigs.py $asm_BLASTN_intron_out $intron_seq_path $asm_path $asm_fix_path
