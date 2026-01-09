#!/bin/bash

for prefix in Train_Ng_0.05 Train_Ng_0.02 Train_Ng_0.01 Train_Ng_0.005
do
    for cov in 30x
    do
        #  Remove intron from assembled contigs
        bash Attached-13-remove_intron_cov.sh $prefix $cov

        #  Remove intron from corrected-reads-assembled contigs
        bash Attached-15-remove_intron_corrected_cov.sh $prefix $cov
    done
done