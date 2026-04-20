#!/bin/bash
prefix=$1
in_fq1_path=$2
in_fq2_path=$3

cov=30x
FQ1_path=fastq/$prefix.$cov.read1.fastq.gz
FQ2_path=fastq/$prefix.$cov.read2.fastq.gz
mkdir fastq 2>/dev/null
mkdir fasta 2>/dev/null
mkdir db 2>/dev/null
#  Sample 1200000 pair-ended reads from the original FastQ files
seqkit sample -n 1200000 -s 100 -o $FQ1_path $in_fq1_path
seqkit sample -n 1200000 -s 100 -o $FQ2_path $in_fq2_path

#  Build BLAST database
for num in 1 2
do
    FQ_path=fastq/$prefix.$cov.read${num}.fastq.gz
    FA_path=fasta/$prefix.$cov.read${num}.fasta
    DB_name=db/$prefix.$cov.db${num}

    #  Convert FASTQ file to FASTA format
    gunzip -c $FQ_path | awk 'NR%4 == 1 {gsub(/^@/, "", $0); printf ">%s\n", $1} NR%4 == 2 {print}' > $FA_path

    #  Build BLAST database from reads in FASTA format
    makeblastdb -in $FA_path -dbtype nucl -out $DB_name
done

#  Align, filter and assemble reads of simple genes (genes without intron)
bash Attached-04-simple_cov_test.sh $prefix $cov

#  Align, filter and assemble reads of genes with intron
bash Attached-05-intron_cov_test.sh $prefix $cov
