#!/bin/bash -l

set -euo pipefail

run_scope=${1:-pilot}
case "$run_scope" in
    pilot)
        array_spec="1-1"
        ;;
    full)
        array_spec="1-61"
        ;;
    *)
        echo "Usage: $0 [pilot|full]" >&2
        exit 1
        ;;
esac

project_path=${JLC_HMM_PROJECT_PATH:-/projects/standard/mfiecas/aron0064/JLC-HMM}
rcode_dir="$project_path/Rcode"
run_dir="$project_path/Routputs"
log_dir="$project_path/LogFiles"
data_dir="$run_dir/Data"
sub_script="$project_path/SubLoopJMHMM.sh"

joint_reference="$data_dir/JMHMMRandInitReducedFitMix5Seed33len96.rda"
two_stage_reference="$data_dir/JMHMMNoSurvRandInitReducedFitMix5Seed33len96.rda"

for required_file in "$joint_reference" "$two_stage_reference" "$sub_script"; do
    if [ ! -f "$required_file" ]; then
        echo "Required NHANES JKn file not found: $required_file" >&2
        exit 1
    fi
done

mkdir -p "$log_dir" "$run_dir"

if ! command -v Rscript >/dev/null 2>&1; then
    module load R/4.4.0-openblas-rocky8
fi


fit_mix_num=5
data_source=nhanes
run_bootstrap=false
init_jitter_scale=0
run_leave_one_out_cv=false
use_hot_start=true
simulation_days=7
num_people=5000
true_mix_num=5
save_reduced_output=true
class_selection_run=true
emission_overlap=low
requested_time=48
requested_mem=20
resampling_method=nhanes_jkn

manifest_file="$run_dir/expected_jobs_nhanes_jkn.tsv"
if [ ! -f "$manifest_file" ]; then
    printf "run_id\tslurm_job_id\tmodel_type\treplicate_id\texpected_file\tarray_spec\trequested_time\trequested_mem\n" > "$manifest_file"
fi
run_id=$(date +"%Y%m%d_%H%M%S")

cd "$project_path"

for model_type in joint two_stage; do
    model_label=Joint
    if [ "$model_type" = "two_stage" ]; then
        model_label=TwoStage
    fi
    job_name="JLC_nhanes_jkn_${model_type}"
    submitted_job_id=$(sbatch --parsable \
        --array="$array_spec" \
        --job-name="$job_name" \
        --time="${requested_time}":00:00 \
        --mem="${requested_mem}"gb \
        --chdir="$project_path" \
        --output="$log_dir/%A_%a.out" \
        --error="$log_dir/%A_%a.err" \
        "$sub_script" \
        "$fit_mix_num" "$model_type" "$data_source" \
        "$run_bootstrap" "$init_jitter_scale" \
        "$run_leave_one_out_cv" "$use_hot_start" \
        "$simulation_days" "$num_people" "$true_mix_num" \
        "$save_reduced_output" "$class_selection_run" \
        "$emission_overlap" "$requested_time" "$requested_mem" \
        "$resampling_method")

    replicate_start=${array_spec%-*}
    replicate_end=${array_spec#*-}
    for replicate_id in $(seq "$replicate_start" "$replicate_end"); do
        expected_file=$(printf \
            "Routputs/JMHMMJKn%sReducedFitMix5RefSeed33Rep%03dlen96.rda" \
            "$model_label" "$replicate_id")
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$run_id" "$submitted_job_id" "$model_type" \
            "$replicate_id" "$expected_file" "$array_spec" \
            "$requested_time" "$requested_mem" >> "$manifest_file"
    done
    echo "Submitted $model_type NHANES JKn array $array_spec as job $submitted_job_id"
done
