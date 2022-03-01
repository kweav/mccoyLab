#!/bin/bash

#SBATCH --job-name=simWL
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --account=rmccoy22

module load R/4.0.2
module load gcc/5.5.0

echo JobID_${SLURM_JOB_ID}_g${1}_s${2}_c${3}_se${4}_ar${5}_r${6}_w${7}_o${8} > /scratch/groups/rmccoy22/kweave23/assess_rhapsodi/vary_window/out_g${1}_s${2}_c${3}_se${4}_ar${5}_r${6}_w${7}_o${8}.txt
date; time Rscript assess_varying_window.R $1 $2 $3 $4 $5 $6 8 $7 $8 &> /scratch/groups/rmccoy22/kweave23/assess_rhapsodi/vary_window/out_g${1}_s${2}_c${3}_se${4}_ar${5}_r${6}_w${7}_o${8}.txt
