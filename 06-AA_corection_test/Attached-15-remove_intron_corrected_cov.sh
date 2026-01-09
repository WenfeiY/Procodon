#!/bin/bash

#  Read prefix and coverage
prefix=$1
cov=$2

#  List of genes with introns
gene_name_list_file=intron_gene_list.txt

#  Default directories of results
root_dir=$prefix

#  Load all gene name
total_commands=$(wc -l < "$protein_name_list_file")

#  20 parallel processing groups
parallel_groups=20

#  Calculate number of tasks for each group
commands_per_group=$(( (total_commands + parallel_groups - 1) / parallel_groups ))

#  Process genes in a parallel multitasking manner
mapfile -t all_protein_name < "$protein_name_list_file"

for ((group=0; group<parallel_groups; group++)); do
    (
        #  Calculate start and end index
        start=$(( group * commands_per_group ))
        end=$(( start + commands_per_group - 1 ))

        #  Run commands
        for ((i=start; i<=end; i++)); do
            if (( i < total_commands )); then
                #  Assign gene name
                protein_name="${all_protein_name[i]}"
                #  Process gene
                echo $protein_name
		        bash Attached-14-remove_intron_corrected_asm.sh $prefix $cov $protein_name
            fi
        done
    ) &
done
