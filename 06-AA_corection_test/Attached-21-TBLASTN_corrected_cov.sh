#!/bin/bash

#  Read prefix and coverage
prefix=$1
cov=$2

#  List of all genes
gene_name_list_file=all_gene_list.txt

#  Default directories of results
root_dir=$prefix
TBLASTN_asm_dir=$root_dir/TBLASTN_asm_out_fix/$cov

#  Create directory
mkdir -p $TBLASTN_asm_dir 2>/dev/null

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
		        bash Attached-19-TBLASTN_corrected_asm.sh $prefix $cov $gene_name
            fi
        done
    ) &
done
