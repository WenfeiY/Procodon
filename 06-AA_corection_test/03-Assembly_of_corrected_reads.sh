#!/bin/bash

for prefix in Train_Ng_0.05 Train_Ng_0.02 Train_Ng_0.01 Train_Ng_0.005
do
    for cov in  30x
    do
        #  Assemble corrected reads
        bash Attached-10-cov_fix_asm.sh $prefix $cov
    done
done
