#!/bin/bash

for ratio in $(seq 0 0.01 0.5)
do
    filename=oHBA_rgs_$ratio.out
    echo "working for $filename"    
    # get SCF ground-state energy
    e_scf=$(grep -i "Total energy in the final basis set =" $filename | awk '{print $9}')
    echo $e_scf > $filename.energy
    grep -A1000 "TDDFT Excitation Energies" oHBA_rgs_$ratio.out | grep -A2 -B2 'Singlet' | grep -i "excitation energy (eV)" | awk '{print $8}' >> $filename.energy
    grep -A1000 "TDDFT Excitation Energies" oHBA_rgs_$ratio.out  | grep -A2 -B2 'Singlet' |  grep -i "Trans. Mom.:" | awk '{print $3 " " $5 " " $7}' > $filename.transdip
done
