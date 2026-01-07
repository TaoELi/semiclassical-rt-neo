#!/bin/bash

#SBATCH --job-name=dtoHBA_pgrad_0_mv_0_rescf_0_dtsm
#SBATCH --time=28-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem -nt 16 oHBA_rgs.in oHBA_rgs.out

