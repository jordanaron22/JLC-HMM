#!/bin/bash -l

TimeArray=("0" "36" "48" "72" "96" "120" "144" "160" "172" "196")
MemArray=("0" "10" "10" "15" "20" "20" "25" "25" "30" "30" )
true_mix_num=5
for model_type in joint two_stage; do
    for fit_mix_num in 1 2 3 4 5 6 7; do
        for data_source in simulation; do
            for run_bootstrap in false; do
                for init_jitter_scale in 1; do
                    for run_leave_one_out_cv in false; do
                        for use_hot_start in false; do
                            for sim_scenario in 0; do
                                sbatch --time="${TimeArray[$fit_mix_num]}":00:00 --mem="${MemArray[$fit_mix_num]}"gb SubLoopJMHMM.sh $fit_mix_num $model_type $data_source $run_bootstrap $init_jitter_scale $run_leave_one_out_cv $use_hot_start $sim_scenario $true_mix_num
                            done
                        done
                    done
                done
            done
        done
    done
done
