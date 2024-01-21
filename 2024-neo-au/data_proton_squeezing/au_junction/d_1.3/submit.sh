#!/bin/bash

#SBATCH --job-name=au
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=32
#SBATCH --partition=general
#SBATCH -p day

qchem -nt 32 au12h.in au12h.out

