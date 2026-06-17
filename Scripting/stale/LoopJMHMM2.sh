#!/bin/bash -l

TimeArray=("0" "36" "48" "72" "96" "36" "144" "160" "172" "196")
MemArray=("0" "10" "10" "15" "20" "20" "25" "25" "30" "30" )
for surv in 2 0; do
    for RE_num in 5; do
        for real_data in 0; do
            for bootstrap in 0; do
                for randomize_init in 0; do
                    for leave_out in 0; do
                        for load_data in 0; do
                            for subset_data in 0; do
                                for sim_size in 0 1 2; do
                                    sbatch --time="${TimeArray[$RE_num]}":00:00 --mem="${MemArray[$RE_num]}"gb SubLoopJMHMM.sh $RE_num $surv $real_data $bootstrap $randomize_init $leave_out $load_data $subset_data $sim_size
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done

subset_data
