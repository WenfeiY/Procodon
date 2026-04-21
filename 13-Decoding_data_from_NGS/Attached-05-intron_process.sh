#!/bin/bash

#  Read prefix and coverage
prefix=$1
in_fq1_path=$2
in_fq2_path=$3

#  List of genes with introns
gene_name_list_file=intron_gene_list.txt

#  Default directories of results
root_dir=$prefix
TBLASTN_dir=$root_dir/TBLASTN_out
BLASTN_dir=$root_dir/BLASTN_out
Read_list_dir=$root_dir/Read_list
FQ_sub_dir=$root_dir/FQ_sub
asm_root_dir=$root_dir/spades.assemble

#  Create directories
mkdir -p $TBLASTN_dir 2>/dev/null
mkdir -p $BLASTN_dir 2>/dev/null
mkdir -p $Read_list_dir 2>/dev/null
mkdir -p $FQ_sub_dir 2>/dev/null
mkdir -p $asm_root_dir 2>/dev/null

#  Load all gene name
total_commands=$(wc -l < "$gene_name_list_file")

#  20 parallel processing groups
parallel_groups=20

#  Calculate number of tasks for each group
commands_per_group=$(( (total_commands + parallel_groups - 1) / parallel_groups ))

#  Process genes in a parallel multitasking manner
mapfile -t all_gene_name < "$gene_name_list_file"

for ((group=0; group<parallel_groups; group++)); do
    (
        #  Calculate start and end index
        start=$(( group * commands_per_group ))
        end=$(( start + commands_per_group - 1 ))

        #  Run commands
        for ((i=start; i<=end; i++)); do
            if (( i < total_commands )); then
                #  Assign gene name
                gene_name="${all_gene_name[i]}"
                #  Process gene
                echo $gene_name
		        bash Attached-03-intron_gene_process.sh $prefix $gene_name $in_fq1_path $in_fq2_path
            fi
        done
    ) &
done



