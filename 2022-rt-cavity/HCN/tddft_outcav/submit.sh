#!/bin/bash

#SBATCH --job-name=HCN_rtdft
#SBATCH --time=21-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem HCN.in HCN.out
