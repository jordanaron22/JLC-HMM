#!/bin/bash -l

declare -A ResourceLimit

array_start=1
array_end=10

# key: simulation_days|num_people|fit_mix_num|model_type|emission_overlap
# value: recommended_time_hours memory_limit_gb
<<<<<<< HEAD
ResourceLimit["1|5000|5|joint|low"]="8 10"
ResourceLimit["3|5000|5|joint|low"]="8 10"
ResourceLimit["7|5000|5|joint|low"]="8 24"
=======
ResourceLimit["1|5000|5|joint|low"]="12 10"
ResourceLimit["3|5000|5|joint|low"]="12 10"
ResourceLimit["7|5000|5|joint|low"]="12 20"
>>>>>>> ae20d82bd6fc9474cd60407d4094256349e9fdbe

ResourceLimit["1|5000|5|joint|high"]="13 10"
ResourceLimit["3|5000|5|joint|high"]="19 10"
ResourceLimit["7|5000|5|joint|high"]="29 20"

<<<<<<< HEAD
ResourceLimit["1|5000|5|two_stage|low"]="8 10"
ResourceLimit["3|5000|5|two_stage|low"]="8 10"
ResourceLimit["7|5000|5|two_stage|low"]="8 24"
=======
ResourceLimit["1|5000|5|two_stage|low"]="12 10"
ResourceLimit["3|5000|5|two_stage|low"]="12 10"
ResourceLimit["7|5000|5|two_stage|low"]="12 20"
>>>>>>> ae20d82bd6fc9474cd60407d4094256349e9fdbe

ResourceLimit["1|5000|5|two_stage|high"]="12 10"
ResourceLimit["3|5000|5|two_stage|high"]="19 10"
ResourceLimit["7|5000|5|two_stage|high"]="28 20"

ResourceLimit["1|5000|2|joint|low"]="8 10"
<<<<<<< HEAD
ResourceLimit["1|5000|3|joint|low"]="8 10"
ResourceLimit["1|5000|4|joint|low"]="8 10"
ResourceLimit["1|5000|6|joint|low"]="10 10"
ResourceLimit["1|5000|7|joint|low"]="10 10"
ResourceLimit["1|5000|8|joint|low"]="10 10"

ResourceLimit["3|5000|2|joint|low"]="8 10"
ResourceLimit["3|5000|3|joint|low"]="8 10"
ResourceLimit["3|5000|4|joint|low"]="8 10"
ResourceLimit["3|5000|6|joint|low"]="16 10"
ResourceLimit["3|5000|7|joint|low"]="18 12"
ResourceLimit["3|5000|8|joint|low"]="22 12"

ResourceLimit["7|5000|2|joint|low"]="10 18"
ResourceLimit["7|5000|3|joint|low"]="12 18"
ResourceLimit["7|5000|4|joint|low"]="17 20"
ResourceLimit["7|5000|6|joint|low"]="35 24"
ResourceLimit["7|5000|7|joint|low"]="38 28"
ResourceLimit["7|5000|8|joint|low"]="37 28"
=======
ResourceLimit["1|5000|3|joint|low"]="15 10"
ResourceLimit["1|5000|4|joint|low"]="15 10"
ResourceLimit["1|5000|6|joint|low"]="20 10"
ResourceLimit["1|5000|7|joint|low"]="20 10"
ResourceLimit["1|5000|8|joint|low"]="20 10"

ResourceLimit["3|5000|2|joint|low"]="15 10"
ResourceLimit["3|5000|3|joint|low"]="15 10"
ResourceLimit["3|5000|4|joint|low"]="15 10"
ResourceLimit["3|5000|6|joint|low"]="20 10"
ResourceLimit["3|5000|7|joint|low"]="30 12"
ResourceLimit["3|5000|8|joint|low"]="30 12"

ResourceLimit["7|5000|2|joint|low"]="30 20"
ResourceLimit["7|5000|3|joint|low"]="30 20"
ResourceLimit["7|5000|4|joint|low"]="30 20"
ResourceLimit["7|5000|6|joint|low"]="45 25"
ResourceLimit["7|5000|7|joint|low"]="45 30"
ResourceLimit["7|5000|8|joint|low"]="50 30"
>>>>>>> ae20d82bd6fc9474cd60407d4094256349e9fdbe

ScenarioKeys=(
    "1|5000|5|joint|low"
    "3|5000|5|joint|low"
    "7|5000|5|joint|low"

    "1|5000|5|joint|high"
    "3|5000|5|joint|high"
    "7|5000|5|joint|high"

    "1|5000|5|two_stage|low"
    "3|5000|5|two_stage|low"
    "7|5000|5|two_stage|low"

    "1|5000|5|two_stage|high"
    "3|5000|5|two_stage|high"
    "7|5000|5|two_stage|high"

    "1|5000|2|joint|low"
    "1|5000|3|joint|low"
    "1|5000|4|joint|low"
    "1|5000|6|joint|low"
    "1|5000|7|joint|low"
    "1|5000|8|joint|low"

    "3|5000|2|joint|low"
    "3|5000|3|joint|low"
    "3|5000|4|joint|low"
    "3|5000|6|joint|low"
    "3|5000|7|joint|low"
    "3|5000|8|joint|low"

    "7|5000|2|joint|low"
    "7|5000|3|joint|low"
    "7|5000|4|joint|low"
    "7|5000|6|joint|low"
    "7|5000|7|joint|low"
    "7|5000|8|joint|low"
)

manifest_file="expected_jobs.tsv"
run_id=$(date +"%Y%m%d_%H%M%S")

mkdir -p "$(dirname "$manifest_file")"

if [ ! -f "$manifest_file" ]; then
    printf "run_id\tslurm_job_id\tjob_name\tsim_num\texpected_file\tfit_mix_num\tmodel_type\tdata_source\trun_bootstrap\tinit_jitter_scale\trun_leave_one_out_cv\tuse_hot_start\tsimulation_days\tnum_people\ttrue_mix_num\tsave_reduced_output\tclass_selection_run\temission_overlap\trequested_time\trequested_mem\n" > "$manifest_file"
fi

true_mix_num=5
save_reduced_output=true
class_selection_run=true
use_hot_start=false
init_jitter_scale=0.1
data_source="simulation"
run_bootstrap=false
run_leave_one_out_cv=false

for scenario_key in "${ScenarioKeys[@]}"; do
    IFS='|' read -r simulation_days num_people fit_mix_num model_type emission_overlap <<< "$scenario_key"

    if [ -z "${ResourceLimit[$scenario_key]+x}" ]; then
        echo "Missing resource limit for scenario: $scenario_key" >&2
        exit 1
    fi

    read -r requested_time requested_mem <<< "${ResourceLimit[$scenario_key]}"

    job_name="JLC_d${simulation_days}_n${num_people}_true${true_mix_num}_fit${fit_mix_num}_${model_type}_${emission_overlap}"

    case "$emission_overlap" in
        low)
            overlap_label="Low"
            ;;
        high)
            overlap_label="High"
            ;;
    esac

    expected_file_prefix="JMHMMDays${simulation_days}People${num_people}Overlap${overlap_label}"

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

    submitted_job_id=$(sbatch --parsable --job-name="$job_name" --time="${requested_time}":00:00 --mem="${requested_mem}"gb SubLoopJMHMM.sh $fit_mix_num $model_type $data_source $run_bootstrap $init_jitter_scale $run_leave_one_out_cv $use_hot_start $simulation_days $num_people $true_mix_num $save_reduced_output $class_selection_run $emission_overlap $requested_time $requested_mem)

    for sim_num in $(seq "$array_start" "$array_end"); do
        expected_file="Routputs/${expected_file_prefix}TrueMix${true_mix_num}FitMix${fit_mix_num}Seed${sim_num}len96.rda"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$run_id" "$submitted_job_id" "$job_name" "$sim_num" "$expected_file" "$fit_mix_num" "$model_type" "$data_source" "$run_bootstrap" "$init_jitter_scale" "$run_leave_one_out_cv" "$use_hot_start" "$simulation_days" "$num_people" "$true_mix_num" "$save_reduced_output" "$class_selection_run" "$emission_overlap" "$requested_time" "$requested_mem" >> "$manifest_file"
    done
done
