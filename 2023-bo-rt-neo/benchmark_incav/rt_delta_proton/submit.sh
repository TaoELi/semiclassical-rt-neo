#!/bin/bash

#SBATCH --job-name=BO_pulse_0.04
#SBATCH --time=14-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem -nt 16 HCN.in HCN.out

