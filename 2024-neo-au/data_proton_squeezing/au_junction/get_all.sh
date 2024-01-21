#!/bin/bash

rm energy.txt dipole.txt
echo "# R(Au-Au)  Escf  Eneo  Eneo(noke)" >> energy.txt
echo "# R(Au-Au)  proton<x>  proton<y>  proton<z>  proton<x^2>  proton<y^2>  proton<z^2>" >> dipole.txt

for dis in 1.4 1.5 1.6 1.62 1.64 1.66 1.665 1.670 1.675 1.68 1.685 1.690 1.695 1.7 1.72 1.74 1.76 1.78 1.8 1.9 2.0 2.1 2.2 2.3 2.4
do


dir=d_$dis

echo "working on $dir"

# 1. get modified protonic densities where the scale is 1000 now
python generate_modified_densities.py $dir
#python generate_diff_densities.py $dir

# 2. get energies

echo "finding energies..."
Escf=$(cat $dir/au12h.out | grep -i "SCF   energy in the final basis set" | awk '{print $9}')
Eneo=$(cat $dir/au12h.out | grep -i "E(NEO-SCF)" | awk '{print $3}')
Enok=$(cat $dir/au12h.out | grep -i "E(NEO-SCF excluding KE quantum protons)" | awk '{print $7}')
echo "$dis  $Escf  $Eneo  $Enok" >> energy.txt

echo "finding dipole moments..."
dipole=$(grep -A1 -i "Protonic dipole moment:" $dir/au12h.out | grep -i "Dx" | awk '{print $2 " " $4 " " $6}')
quadrupole=$(grep -A1 -i "Protonic quadrupole moment:" $dir/au12h.out | grep -i "Qx" | awk '{print $2 " " $4 " " $6}')
echo "$dis $dipole $quadrupole" >> dipole.txt

done

