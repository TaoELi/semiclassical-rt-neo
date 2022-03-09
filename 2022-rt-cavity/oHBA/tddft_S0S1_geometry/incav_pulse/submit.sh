#!/bin/bash

#SBATCH --job-name=HCN_rtdft
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH -p day

qchem -nt 16 oHBA_rgs.in oHBA_rgs.out
