#!/bin/bash

#SBATCH --job-name=HCN_rtdft
#SBATCH --time=21-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem -nt 4 HCN.in HCN.out
