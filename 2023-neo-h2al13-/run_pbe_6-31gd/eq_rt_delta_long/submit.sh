#!/bin/bash

#SBATCH --job-name=Al13_cl_mv
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=32
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem -nt 32 Al13H2.in Al13H2.out

