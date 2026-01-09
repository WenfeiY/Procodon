#!/bin/bash

#  Read prefix and coverage
prefix=$1
cov=$2

#  List of simple genes
gene_name_list_file=simple_gene_list.txt

#  Default directories of results
root_dir=$prefix
TBLASTN_dir=$root_dir/TBLASTN_out/$cov
BLASTN_dir=$root_dir/BLASTN_out/$cov
Read_list_dir=$root_dir/Read_list/$cov
FQ_sub_dir=$root_dir/FQ_sub/$cov
asm_root_dir=$root_dir/spades.assemble/$cov

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
		        bash Attached-02-simple_gene_cov_process.sh $prefix $cov $gene_name
            fi
        done
    ) &
done
