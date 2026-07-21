#!/bin/bash -l

declare -A ResourceLimit

full_run_array_start=11
full_run_array_end=110
cv_fold_count=20

# key: fit_mix_num|model_type
# value: recommended_time_hours memory_limit_gb

ResourceLimit["1|joint"]="6 5"
ResourceLimit["2|joint"]="12 6"
ResourceLimit["3|joint"]="12 7"
ResourceLimit["4|joint"]="12 8"
ResourceLimit["5|joint"]="12 9"
ResourceLimit["6|joint"]="12 10"
ResourceLimit["7|joint"]="12 11"
ResourceLimit["8|joint"]="12 12"
ResourceLimit["9|joint"]="12 14"

ResourceLimit["1|two_stage"]="6 5"
ResourceLimit["2|two_stage"]="12 6"
ResourceLimit["3|two_stage"]="12 7"
ResourceLimit["4|two_stage"]="12 8"
ResourceLimit["5|two_stage"]="12 9"
ResourceLimit["6|two_stage"]="12 10"
ResourceLimit["7|two_stage"]="12 11"
ResourceLimit["8|two_stage"]="12 12"
ResourceLimit["9|two_stage"]="12 14"


ScenarioKeys=(
    "1|joint"
    # "2|joint"
    # "3|joint"
    # "4|joint"
    "5|joint"
    # "6|joint"
    # "7|joint"
    # "8|joint"
    # "9|joint"

    "1|two_stage"
    # "2|two_stage"
    # "3|two_stage"
    # "4|two_stage"
    "5|two_stage"
    # "6|two_stage"
    # "7|two_stage"
    # "8|two_stage"
    # "9|two_stage"
)

manifest_file="expected_jobs_nhanes.tsv"
run_id=$(date +"%Y%m%d_%H%M%S")

mkdir -p "$(dirname "$manifest_file")"

if [ ! -f "$manifest_file" ]; then
    printf "run_id\tslurm_job_id\tjob_name\tsim_num\texpected_file\tfit_mix_num\tmodel_type\tdata_source\trun_bootstrap\tinit_jitter_scale\trun_leave_one_out_cv\tuse_hot_start\tsimulation_days\tnum_people\ttrue_mix_num\tsave_reduced_output\tclass_selection_run\temission_overlap\trequested_time\trequested_mem\n" > "$manifest_file"
fi

data_source="nhanes"
run_bootstrap=false
init_jitter_scale=0.05
run_leave_one_out_cv=true
use_hot_start=true
simulation_days=7
num_people=5000
save_reduced_output=false
class_selection_run=false
emission_overlap="low"

if [ "$run_leave_one_out_cv" = "true" ]; then
    # For NHANES CV, SLURM_ARRAY_TASK_ID is the held-out fold_id.
    array_start=1
    array_end=$cv_fold_count
else
    array_start=$full_run_array_start
    array_end=$full_run_array_end
fi

for scenario_key in "${ScenarioKeys[@]}"; do
    IFS='|' read -r fit_mix_num model_type <<< "$scenario_key"

    if [ -z "${ResourceLimit[$scenario_key]+x}" ]; then
        echo "Missing resource limit for scenario: $scenario_key" >&2
        exit 1
    fi

    read -r requested_time requested_mem <<< "${ResourceLimit[$scenario_key]}"

    true_mix_num="$fit_mix_num"
    job_name="JLC_nhanes_fit${fit_mix_num}_${model_type}"
    expected_file_prefix="JMHMM"

    if [ "$model_type" = "two_stage" ]; then
        expected_file_prefix="${expected_file_prefix}NoSurv"
    fi

    if [ "$run_bootstrap" = "true" ]; then
        expected_file_prefix="${expected_file_prefix}Bootstrap"
    fi

    if [ "$run_leave_one_out_cv" = "true" ]; then
        expected_file_prefix="${expected_file_prefix}LeaveOut"
    fi

    if [ "$init_jitter_scale" != "0" ] && [ "$init_jitter_scale" != "0.0" ]; then
        expected_file_prefix="${expected_file_prefix}RandInit"
    fi

    if [ "$use_hot_start" = "true" ]; then
        expected_file_prefix="${expected_file_prefix}LoadIn"
    fi

    if [ "$class_selection_run" = "true" ]; then
        expected_file_prefix="${expected_file_prefix}ClassSelection"
    fi

    if [ "$save_reduced_output" = "true" ]; then
        expected_file_prefix="${expected_file_prefix}Reduced"
    fi

    submitted_job_id=$(sbatch --parsable --array="${array_start}-${array_end}" --job-name="$job_name" --time="${requested_time}":00:00 --mem="${requested_mem}"gb SubLoopJMHMM.sh $fit_mix_num $model_type $data_source $run_bootstrap $init_jitter_scale $run_leave_one_out_cv $use_hot_start $simulation_days $num_people $true_mix_num $save_reduced_output $class_selection_run $emission_overlap $requested_time $requested_mem)

    for sim_num in $(seq "$array_start" "$array_end"); do
        expected_file="Routputs/${expected_file_prefix}FitMix${fit_mix_num}Seed${sim_num}len96.rda"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$run_id" "$submitted_job_id" "$job_name" "$sim_num" "$expected_file" "$fit_mix_num" "$model_type" "$data_source" "$run_bootstrap" "$init_jitter_scale" "$run_leave_one_out_cv" "$use_hot_start" "$simulation_days" "$num_people" "$true_mix_num" "$save_reduced_output" "$class_selection_run" "$emission_overlap" "$requested_time" "$requested_mem" >> "$manifest_file"
    done
done
