#!/bin/bash
prefix=$1
in_fq1_path=$2
in_fq2_path=$3

FA1_path=fasta/$prefix.read1.fasta
FA2_path=fasta/$prefix.read2.fasta
DB1_name=db/$prefix.db1
DB2_name=db/$prefix.db2

mkdir fasta 2>/dev/null
mkdir db 2>/dev/null


#  Build BLAST database

#  Convert FASTQ file to FASTA format
gunzip -c $in_fq1_path | awk 'NR%4 == 1 {gsub(/^@/, "", $0); printf ">%s\n", $1} NR%4 == 2 {print}' > $FA1_path
gunzip -c $in_fq2_path | awk 'NR%4 == 1 {gsub(/^@/, "", $0); printf ">%s\n", $1} NR%4 == 2 {print}' > $FA2_path

#  Build BLAST database from reads in FASTA format
makeblastdb -in $FA1_path -dbtype nucl -out $DB1_name
makeblastdb -in $FA2_path -dbtype nucl -out $DB2_name

#  Align, filter and assemble reads of simple genes (genes without intron)
bash Attached-04-simple_process.sh $prefix $in_fq1_path $in_fq2_path

#  Align, filter and assemble reads of genes with intron
bash Attached-05-intron_process.sh $prefix $in_fq1_path $in_fq2_path
