#!/bin/bash

#  Read prefix, coverage and gene name
prefix=$1
cov=$2
gene=$3

#  The directory of protein sequences from yeast genome
protein_root_dir="./S.cerevisiae.Protein"
#  The directory of gene flanking sequences (upstream and downstream 50 bp) from yeast genome
flanking_root_dir="./S.cerevisiae.Gene_flanking"

#  Default directories of results
root_dir=$prefix
TBLASTN_dir=$root_dir/TBLASTN_out/$cov
BLASTN_dir=$root_dir/BLASTN_out/$cov
Read_list_dir=$root_dir/Read_list/$cov
FQ_sub_dir=$root_dir/FQ_sub/$cov
asm_root_dir=$root_dir/spades.assemble/$cov

#  Input NGS data path
FQ1_path=fastq/$prefix.$cov.read1.fastq.gz
FQ2_path=fastq/$prefix.$cov.read2.fastq.gz
FA1_path=fasta/$prefix.$cov.read1.fasta
FA2_path=fasta/$prefix.$cov.read2.fasta
DB1_name=db/$prefix.$cov.db1
DB2_name=db/$prefix.$cov.db2

#  FASTA file path of protein sequence
protein_seq_path=$protein_root_dir/$gene.Protein.fa
#  FASTA file path of flanking sequences
flanking_seq_path=$flanking_root_dir/$gene.flanking.fa

#  TBLASTN results of reads aligned to protein sequences
TBLASTN1_out=$TBLASTN_dir/$gene.tblastn.1.out
TBLASTN2_out=$TBLASTN_dir/$gene.tblastn.2.out
#  BLASTN results of reads aligned to flanking sequences
BLASTN1_out=$BLASTN_dir/$gene.blastn.1.out
BLASTN2_out=$BLASTN_dir/$gene.blastn.2.out

#  TBLASTN (protein sequences to pair-ended reads)
tblastn -query $protein_seq_path \
    -db $DB1_name \
    -outfmt 6 \
    -evalue 1e-5 \
    -matrix BLOSUM80 \
    -soft_masking true \
    -max_target_seqs 500000 \
    > $TBLASTN1_out

tblastn -query $protein_seq_path \
    -db $DB2_name \
    -outfmt 6 \
    -evalue 1e-5 \
    -matrix BLOSUM80 \
    -soft_masking true \
    -max_target_seqs 500000 \
    > $TBLASTN2_out

#  File path of Filtered TBLASTN results
fil_TBLASTN1_out=$TBLASTN_dir/$gene.tblastn.filt.1.out
fil_TBLASTN2_out=$TBLASTN_dir/$gene.tblastn.filt.2.out

#  Filter TBLASTN results by bit-score
python Attached-01-filter_TBLASTN.py $protein_seq_path $TBLASTN1_out $fil_TBLASTN1_out
python Attached-01-filter_TBLASTN.py $protein_seq_path $TBLASTN2_out $fil_TBLASTN2_out

#BLASTN align flanking regions and reads
blastn -query $flanking_seq_path \
    -db $DB1_name \
    -outfmt 6 \
    -word_size 15 \
    | awk '$12 > 65' \
    > $BLASTN1_out

blastn -query $flanking_seq_path \
    -db $DB2_name \
    -outfmt 6 \
    -word_size 15 \
    | awk '$12 > 65' \
    > $BLASTN2_out

#  File path of reads list
Read_list=$Read_list_dir/$gene.read.list.txt
Read_list1=$Read_list_dir/$gene.read.list1.txt
Read_list2=$Read_list_dir/$gene.read.list2.txt

#  Extract reads ID, remove duplicates
cut -f2 $fil_TBLASTN1_out | sort -u > $Read_list && \
    cut -f2 $fil_TBLASTN2_out | sort -u >> $Read_list && \
    cut -f2 $BLASTN1_out | sort -u >> $Read_list && \
    cut -f2 $BLASTN2_out | sort -u >> $Read_list && \
    sort -u $Read_list -o $Read_list
wait

#  Add surfix
awk '{print substr($0, 1, length($0)-2)}' "$Read_list" | sort -u | tee \
    >(awk -v outfile="$Read_list1" '{print $0 "/1" > outfile}') \
    >(awk -v outfile="$Read_list2" '{print $0 "/2" > outfile}') >/dev/null
wait

#  File path of FASTQ subset
FQ1_sub=$FQ_sub_dir/$gene.read1.fastq
FQ2_sub=$FQ_sub_dir/$gene.read2.fastq

#  Extract FASTQ subset
seqkit grep -f $Read_list1 $FQ1_path -o $FQ1_sub
seqkit grep -f $Read_list2 $FQ2_path -o $FQ2_sub

#  Assemble reads
spades_dir=$asm_root_dir/$gene
spades.py -o $spades_dir/ -1 $FQ1_sub -2 $FQ2_sub --careful -k 21,33,55 -t 8 
