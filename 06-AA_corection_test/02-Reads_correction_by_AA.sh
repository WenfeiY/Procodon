#!/bin/bash

for prefix in Train_Ng_0.05 Train_Ng_0.02 Train_Ng_0.01 Train_Ng_0.005
do
    for cov in 30x
    do
        #  Correct reads by AA
        bash Attached-08-NGS_cov_fix.sh $prefix $cov
    done
done
