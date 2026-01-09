#!/bin/bash

#  Generate 1200000 pair-ended reads for each mutation rate
dwgsim -N 1200000 -1 150 -2 150 -y 0 -e 0.05  S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.fasta Train_Ng_0.05
dwgsim -N 1200000 -1 150 -2 150 -y 0 -e 0.02  S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.fasta Train_Ng_0.02
dwgsim -N 1200000 -1 150 -2 150 -y 0 -e 0.01  S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.fasta Train_Ng_0.01
dwgsim -N 1200000 -1 150 -2 150 -y 0 -e 0.005  S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.fasta Train_Ng_0.005

mkdir fastq 2>/dev/null
mkdir fasta 2>/dev/null
mkdir db 2>/dev/null

#  Change file names and move to directory ./fastq/
for prefix in Train_Ng_0.05 Train_Ng_0.02 Train_Ng_0.01 Train_Ng_0.005
do
    for num in 1 2
    do
        mv $prefix.bwa.read${num}.fastq.gz fastq/$prefix.30x.read${num}.fastq.gz
    done
done

#  Build BLAST database
for prefix in Train_Ng_0.05 Train_Ng_0.02 Train_Ng_0.01 Train_Ng_0.005
do
    for num in 1 2
    do
        for cov in 30x
        do
            FQ_path=fastq/$prefix.$cov.read${num}.fastq.gz
            FA_path=fasta/$prefix.$cov.read${num}.fasta
            DB_name=db/$prefix.$cov.db${num}

            #  Convert FASTQ file to FASTA format
            gunzip -c $FQ_path | awk 'NR%4 == 1 {gsub(/^@/, "", $0); printf ">%s\n", $1} NR%4 == 2 {print}' > $FA_path

            #  Build BLAST database from reads in FASTA format
            makeblastdb -in $FA_path -dbtype nucl -out $DB_name
        done
    done
done

for prefix in Train_Ng_0.05 Train_Ng_0.02 Train_Ng_0.01 Train_Ng_0.005
do
    for cov in 30x
    do
        
        #  Align, filter and assemble reads of simple genes (genes without intron)
        bash Attached-04-simple_cov_test.sh $prefix $cov

        #  Align, filter and assemble reads of genes with intron
        bash Attached-05-intron_cov_test.sh $prefix $cov
    done
done
