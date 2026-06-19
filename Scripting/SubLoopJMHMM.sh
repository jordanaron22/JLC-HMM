#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=NONE
#SBATCH --array=1-10
#SBATCH -o LogFiles/%A_%a.out
#SBATCH -e LogFiles/%A_%a.err
fit_mix_num=$1
model_type=$2
data_source=$3
run_bootstrap=$4
init_jitter_scale=$5
run_leave_one_out_cv=$6
use_hot_start=$7
sim_scenario=$8
true_mix_num=${9:-$fit_mix_num}
save_reduced_output=${10:-false}
class_selection_run=${11:-false}
emission_overlap=${12:-low}
date
path="$HOME/JM"
cd $path/Routputs
start_time=$(date +%s)
module load R/4.4.0
Rscript $path/Rcode/JMHMM.R $fit_mix_num $model_type $data_source $run_bootstrap $init_jitter_scale $run_leave_one_out_cv $use_hot_start $sim_scenario $true_mix_num $save_reduced_output $class_selection_run $emission_overlap
finish_time=$(date +%s)
elapsed_time=$((finish_time  - start_time))

((sec=elapsed_time%60, elapsed_time/=60, min=elapsed_time%60, hrs=elapsed_time/60))
timestamp=$(printf "Total time taken - %d hours, %d minutes, and %d seconds." $hrs $min $sec)
echo $timestamp
