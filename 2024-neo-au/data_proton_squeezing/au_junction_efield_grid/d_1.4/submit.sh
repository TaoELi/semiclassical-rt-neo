#!/bin/bash

#SBATCH --job-name=au
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=32
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem -nt 32 au12h.in au12h.out

