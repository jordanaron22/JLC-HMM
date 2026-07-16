#!/bin/bash -l

declare -A ResourceLimit

array_start=1
array_end=200

# key: simulation_days|num_people|fit_mix_num|model_type|emission_overlap
# value: recommended_time_hours memory_limit_gb
ResourceLimit["1|5000|5|joint|low"]="4 4"
ResourceLimit["3|5000|5|joint|low"]="16 10"
ResourceLimit["7|5000|5|joint|low"]="30 20"

ResourceLimit["1|5000|5|joint|high"]="5 5"
ResourceLimit["3|5000|5|joint|high"]="16 10"
ResourceLimit["7|5000|5|joint|high"]="52 20"

ResourceLimit["1|5000|5|joint|mid"]="5 5"
ResourceLimit["3|5000|5|joint|mid"]="16 10"
ResourceLimit["7|5000|5|joint|mid"]="52 20"


ScenarioKeys=(
    "1|5000|5|joint|low|profiled"
    "3|5000|5|joint|low|profiled"
    "7|5000|5|joint|low|profiled"

    "1|5000|5|joint|mid|profiled"
    "3|5000|5|joint|mid|profiled"
    "7|5000|5|joint|mid|profiled"

    "1|5000|5|joint|high|profiled"
    "3|5000|5|joint|high|profiled"
    "7|5000|5|joint|high|profiled"
)

manifest_file="expected_oakes_jobs.tsv"
run_id=$(date +"%Y%m%d_%H%M%S")
h2_eps="1e-5"

mkdir -p "$(dirname "$manifest_file")"

if [ ! -f "$manifest_file" ]; then
    printf "run_id\tslurm_job_id\tjob_name\tsim_num\tinput_file\texpected_file\tfit_mix_num\tmodel_type\tdata_source\trun_bootstrap\tinit_jitter_scale\trun_leave_one_out_cv\tuse_hot_start\tsimulation_days\tnum_people\ttrue_mix_num\tsave_reduced_output\tclass_selection_run\temission_overlap\tsurvival_baseline_mode\trequested_time\trequested_mem\th2_eps\n" > "$manifest_file"
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
    IFS='|' read -r simulation_days num_people fit_mix_num model_type emission_overlap survival_baseline_mode <<< "$scenario_key"
    survival_baseline_mode=${survival_baseline_mode:-profiled}

    if [ "$model_type" != "joint" ]; then
        echo "Skipping non-joint scenario: $scenario_key"
        continue
    fi

    case "$survival_baseline_mode" in
        fixed|profiled)
            ;;
        *)
            echo "Unknown survival_baseline_mode: $survival_baseline_mode" >&2
            exit 1
            ;;
    esac

    resource_key="${simulation_days}|${num_people}|${fit_mix_num}|${model_type}|${emission_overlap}"
    if [ -z "${ResourceLimit[$resource_key]+x}" ]; then
        echo "Missing resource limit for scenario: $resource_key" >&2
        exit 1
    fi

    read -r requested_time requested_mem <<< "${ResourceLimit[$resource_key]}"

    job_name="JLCse_d${simulation_days}_n${num_people}_true${true_mix_num}_fit${fit_mix_num}_${model_type}_${emission_overlap}_${survival_baseline_mode}"

    case "$emission_overlap" in
        low)
            overlap_label="Low"
            ;;
        mid)
            overlap_label="Mid"
            ;;
        high)
            overlap_label="High"
            ;;
        *)
            echo "Unknown emission_overlap: $emission_overlap" >&2
            exit 1
            ;;
    esac

    expected_file_prefix="JMHMMDays${simulation_days}People${num_people}Overlap${overlap_label}"

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

    submitted_job_id=$(sbatch --parsable --array="${array_start}-${array_end}" --job-name="$job_name" --time="${requested_time}":00:00 --mem="${requested_mem}"gb SubLoopJLC-HMMse.sh $fit_mix_num $model_type $data_source $run_bootstrap $init_jitter_scale $run_leave_one_out_cv $use_hot_start $simulation_days $num_people $true_mix_num $save_reduced_output $class_selection_run $emission_overlap $requested_time $requested_mem $survival_baseline_mode)

    for sim_num in $(seq "$array_start" "$array_end"); do
        input_file="Routputs/${expected_file_prefix}TrueMix${true_mix_num}FitMix${fit_mix_num}Seed${sim_num}len96.rda"
        expected_file="Routputs/${expected_file_prefix}TrueMix${true_mix_num}FitMix${fit_mix_num}Seed${sim_num}len96_oakes_${survival_baseline_mode}_se.rda"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$run_id" "$submitted_job_id" "$job_name" "$sim_num" "$input_file" "$expected_file" "$fit_mix_num" "$model_type" "$data_source" "$run_bootstrap" "$init_jitter_scale" "$run_leave_one_out_cv" "$use_hot_start" "$simulation_days" "$num_people" "$true_mix_num" "$save_reduced_output" "$class_selection_run" "$emission_overlap" "$survival_baseline_mode" "$requested_time" "$requested_mem" "$h2_eps" >> "$manifest_file"
    done
done
