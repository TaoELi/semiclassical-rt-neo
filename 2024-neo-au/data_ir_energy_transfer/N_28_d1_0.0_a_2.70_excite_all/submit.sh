#!/bin/bash

#SBATCH --job-name=au_spacing
#SBATCH --time=3-08:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=32
#SBATCH --partition=general
#SBATCH -p pi_hammes_schiffer

qchem -nt 32 hfau2n.in hfau2n.out

