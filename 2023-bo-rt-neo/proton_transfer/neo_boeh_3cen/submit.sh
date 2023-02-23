#!/bin/bash

#SBATCH --job-name=3cen_boeh_plotting
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem -nt 16 mda.in mda.out

