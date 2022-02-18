#!/bin/bash

# ./run_wgets.sh nl_sumstat_wgets.txt

#SBATCH --job-name=download_sumstats
#SBATCH --time=2:00:00
#SBATCH -N 1

cat $1 | while read line
do
   $line
done
