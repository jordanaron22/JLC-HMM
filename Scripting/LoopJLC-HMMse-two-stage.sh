#!/bin/bash -l

declare -A ResourceLimit

array_start=100
array_end=100

# key: simulation_days|num_people|fit_mix_num|model_type|emission_overlap
# value: recommended_time_hours memory_limit_gb
# These initially mirror the joint-SE requests. Adjust them after reviewing
# the first completed Murphy-Topel jobs for each scenario.
ResourceLimit["1|5000|5|two_stage|low"]="4 4"
ResourceLimit["3|5000|5|two_stage|low"]="16 10"
ResourceLimit["7|5000|5|two_stage|low"]="30 20"

ResourceLimit["1|5000|5|two_stage|high"]="5 5"
ResourceLimit["3|5000|5|two_stage|high"]="16 10"
ResourceLimit["7|5000|5|two_stage|high"]="52 20"

ResourceLimit["1|5000|5|two_stage|low"]="1 4"

ScenarioKeys=(
    "1|5000|5|two_stage|low"
    # "3|5000|5|two_stage|low"
    # "7|5000|5|two_stage|low"

    # "1|5000|5|two_stage|high"
    # "3|5000|5|two_stage|high"
    # "7|5000|5|two_stage|high"
)

manifest_file="expected_two_stage_mt_jobs.tsv"
run_id=$(date +"%Y%m%d_%H%M%S")
h2_eps="1e-5"
# format_eps_label(1e-5) in JLC-HMMse-two-stage.R produces 1e-05.
h2_eps_file_label="1e-05"
retain_score_matrices=false

mkdir -p "$(dirname "$manifest_file")"

if [ ! -f "$manifest_file" ]; then
    printf "run_id\tslurm_job_id\tjob_name\tsim_num\tinput_file\texpected_file\tfit_mix_num\tmodel_type\tdata_source\trun_bootstrap\tinit_jitter_scale\trun_leave_one_out_cv\tuse_hot_start\tsimulation_days\tnum_people\ttrue_mix_num\tsave_reduced_output\tclass_selection_run\temission_overlap\trequested_time\trequested_mem\th2_eps\tretain_score_matrices\n" > "$manifest_file"
fi

#these need to match the load-in files, as they determine file name
#doesnt jitter before se calculation
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

    if [ "$model_type" != "two_stage" ]; then
        echo "Skipping non-two-stage scenario: $scenario_key"
        continue
    fi

    resource_key="${simulation_days}|${num_people}|${fit_mix_num}|${model_type}|${emission_overlap}"
    if [ -z "${ResourceLimit[$resource_key]+x}" ]; then
        echo "Missing resource limit for scenario: $resource_key" >&2
        exit 1
    fi

    read -r requested_time requested_mem <<< "${ResourceLimit[$resource_key]}"

    job_name="JLCmt_d${simulation_days}_n${num_people}_true${true_mix_num}_fit${fit_mix_num}_${emission_overlap}_eps${h2_eps_file_label}"

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

    # build_model_name() places NoSurv before RandInit/LoadIn/ClassSelection.
    expected_file_prefix="${expected_file_prefix}NoSurv"

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

    submitted_job_id=$(sbatch --parsable --array="${array_start}-${array_end}" --job-name="$job_name" --time="${requested_time}":00:00 --mem="${requested_mem}"gb SubLoopJLC-HMMse-two-stage.sh "$fit_mix_num" "$model_type" "$data_source" "$run_bootstrap" "$init_jitter_scale" "$run_leave_one_out_cv" "$use_hot_start" "$simulation_days" "$num_people" "$true_mix_num" "$save_reduced_output" "$class_selection_run" "$emission_overlap" "$requested_time" "$requested_mem" "$h2_eps" "$retain_score_matrices")

    for sim_num in $(seq "$array_start" "$array_end"); do
        input_file="Routputs/${expected_file_prefix}TrueMix${true_mix_num}FitMix${fit_mix_num}Seed${sim_num}len96.rda"
        expected_file="Routputs/${expected_file_prefix}TrueMix${true_mix_num}FitMix${fit_mix_num}Seed${sim_num}len96_two_stage_murphy_topel_eps${h2_eps_file_label}_se.rda"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$run_id" "$submitted_job_id" "$job_name" "$sim_num" "$input_file" "$expected_file" "$fit_mix_num" "$model_type" "$data_source" "$run_bootstrap" "$init_jitter_scale" "$run_leave_one_out_cv" "$use_hot_start" "$simulation_days" "$num_people" "$true_mix_num" "$save_reduced_output" "$class_selection_run" "$emission_overlap" "$requested_time" "$requested_mem" "$h2_eps" "$retain_score_matrices" >> "$manifest_file"
    done
done
