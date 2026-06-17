#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=NONE
#SBATCH --array=1-100
#SBATCH -o LogFiles/%A_%a.out
#SBATCH -e LogFiles/%A_%a.err
RE_num=$1
surv=$2
real_data=$3
bootstrap=$4
randomize_init=$5
leave_out=$6
load_data=$7
subset_data=$8
sim_size=$9
date
path="$HOME/JM"
cd $path/Routputs
start_time=$(date +%s)
module load R/4.4.0
Rscript $path/Rcode/JMHMM.R $RE_num $surv $real_data $bootstrap $randomize_init $leave_out $load_data $subset_data $sim_size
finish_time=$(date +%s)
elapsed_time=$((finish_time  - start_time))

((sec=elapsed_time%60, elapsed_time/=60, min=elapsed_time%60, hrs=elapsed_time/60))
timestamp=$(printf "Total time taken - %d hours, %d minutes, and %d seconds." $hrs $min $sec)
echo $timestamp

