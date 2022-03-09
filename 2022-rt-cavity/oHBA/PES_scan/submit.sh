#!/bin/bash

#SBATCH --job-name=oHBA_LR
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem oHBA_rgs.in oHBA_rgs.out
