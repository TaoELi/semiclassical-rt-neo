#!/bin/bash

for Freq in 0.3475
do
for Loss in 6e4 3e4 8e3 2e3
do
for Coupling in 4e-4
do

	echo "generating a QChem submission file with Freq = $Freq, cavity coupling = $Coupling, and cavity loss lifetime = $Loss"
	dir=Freq_$Freq\_CavLifetime_$Loss\_Coupling_$Coupling
	mkdir $dir	
	cp HCN.in submit.sh $dir

	cd $dir
	
	# modify the input parameters
	sed -i "s/in_cavity true 2.01 2 1e-2 100/in_cavity true $Freq 0 $Coupling $Loss/" HCN.in
	sed -i "s/HCN_rtdft/HCN_rtdft_$dir/" submit.sh
	
	sbatch submit.sh
	cd ..

done    
done
done
