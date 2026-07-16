#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=aron0064@umn.edu
#SBATCH --array=1-200
#SBATCH -A mfiecas
#SBATCH -o LogFiles/%A_%a.out
#SBATCH -e LogFiles/%A_%a.err
fit_mix_num=$1
model_type=$2
data_source=$3
run_bootstrap=$4
init_jitter_scale=$5
run_leave_one_out_cv=$6
use_hot_start=$7
simulation_days=$8
num_people=$9
true_mix_num=${10:-$fit_mix_num}
save_reduced_output=${11:-false}
class_selection_run=${12:-false}
emission_overlap=${13:-low}
time_limit_hours=${14:-NA}
memory_limit_gb=${15:-NA}
survival_baseline_mode=${16:-profiled}

date
path="/projects/standard/mfiecas/aron0064/JLC-HMM"
cd $path/Routputs
start_time=$(date +%s)
module load R/4.4.0-openblas-rocky8
echo "JLC-HMMse survival baseline mode: $survival_baseline_mode"
Rscript $path/Rcode/JLC-HMMse.R $fit_mix_num $model_type $data_source $run_bootstrap $init_jitter_scale $run_leave_one_out_cv $use_hot_start $simulation_days $num_people $true_mix_num $save_reduced_output $class_selection_run $emission_overlap $time_limit_hours $memory_limit_gb $survival_baseline_mode
finish_time=$(date +%s)
elapsed_time=$((finish_time  - start_time))

((sec=elapsed_time%60, elapsed_time/=60, min=elapsed_time%60, hrs=elapsed_time/60))
timestamp=$(printf "Total time taken - %d hours, %d minutes, and %d seconds." $hrs $min $sec)
echo $timestamp
