#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=NONE
#SBATCH --array=1-73
#SBATCH -o LogFiles/%A_%a.out
#SBATCH -e LogFiles/%A_%a.err
RE_num=$1
path="$HOME/JM"
cd $path/Routputs
start_time=$(date +%s)
module load R/4.4.0
Rscript $path/Rcode/JMHMMsim.R $RE_num
finish_time=$(date +%s)
elapsed_time=$((finish_time  - start_time))

((sec=elapsed_time%60, elapsed_time/=60, min=elapsed_time%60, hrs=elapsed_time/60))
timestamp=$(printf "Total time taken - %d hours, %d minutes, and %d seconds." $hrs $min $sec)
echo $timestamp

