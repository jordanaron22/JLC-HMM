#!/bin/bash -l

TimeArray=("1" "1" "10" "24")
MemArray=("5" "5" "10" "15")
for RE_num in 3; do
    sbatch --time="${TimeArray[$RE_num - 1]}":00:00 --mem="${MemArray[$RE_num - 1]}"gb JMHMMscript.sh $RE_num 
done

