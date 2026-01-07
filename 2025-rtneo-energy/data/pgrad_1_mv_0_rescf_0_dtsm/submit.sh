#!/bin/bash

#SBATCH --job-name=oHBA_pgrad_0_mv_0_rescf_0
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH -p week

qchem -nt 16 oHBA_rgs.in oHBA_rgs.out

