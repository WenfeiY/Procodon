#!/bin/bash

for prefix in Train_Ng_0.05 Train_Ng_0.02 Train_Ng_0.01 Train_Ng_0.005
do
    for cov in  30x
    do
        #  Align contigs (without correction) with protein sequences
        bash TBLASTN_cov.sh $prefix $cov

        #  Align contigs (corrected) with protein sequences
        bash TBLASTN_cov_fix.sh $prefix $cov
    done
done
