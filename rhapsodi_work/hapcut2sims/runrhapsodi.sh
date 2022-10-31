#!/bin/bash

#SBATCH --partition=defq
#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=185GB
#SBATCH --account=rmccoy22
#SBATCH --array=2918-5832

ml r/4.0.2

arg_arr=($(sed "${SLURM_ARRAY_TASK_ID}q;d" conditions.txt))
#arg_arr=($(sed "2917q;d" conditions.txt)) #used for testing

# ${arg_arr[0]} is number of gametes
# ${arg_arr[1]} is number of snps
# ${arg_arr[2]} is coverage
# ${arg_arr[3]} is seq error rate
# ${arg_arr[4]} is avg recomb rate
# ${arg_arr[5]} is random seed

filelocation_base1a="/home/kweave23/data_rmccoy22/sperm_seq_rhapsodi/genModel_output/mse_0.005_mar_1/"
filelocation_base2a="g${arg_arr[0]}_s${arg_arr[1]}_c${arg_arr[2]}_se${arg_arr[3]}_r${arg_arr[4]}"

filelocation_base1b="/home/kweave23/data_rmccoy22/sperm_seq_rhapsodi/genModel_output/changing_model_params/"
filelocation_base2b="g${arg_arr[0]}_s${arg_arr[1]}_c${arg_arr[2]}_se${arg_arr[3]}_r${arg_arr[4]}_mse${arg_arr[3]}_mr${arg_arr[4]}"

genfile_base="runGen_gam_${arg_arr[0]}_snp_${arg_arr[1]}_cov_${arg_arr[2]}_seqerr_${arg_arr[3]}_avgr_${arg_arr[4]}_rs_${arg_arr[5]}"
sim_base="${filelocation_base1a}${filelocation_base2a}/${genfile_base}"
sim_out_dir="${filelocation_base1b}${filelocation_base2b}"


if ( [ ${arg_arr[3]} != "0.005" ] && [ ${arg_arr[4]} != "1" ] ); then
  Rscript assess_with_rhapsodi_mparams.R ${arg_arr[0]} ${arg_arr[1]} ${arg_arr[2]} ${arg_arr[3]} ${arg_arr[4]} ${arg_arr[5]} 48 $sim_base $sim_out_dir
fi
