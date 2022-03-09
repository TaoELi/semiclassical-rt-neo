#!/bin/bash

for Freq in 13.3344764488
do
for Loss in 1e2
do
for Coupling in 1e-3 2e-3 3e-3 5e-3 6e-3 
do

	echo "generating a QChem submission file with Freq = $Freq, cavity coupling = $Coupling, and cavity loss lifetime = $Loss"
	dir=Freq_$Freq\_CavLifetime_$Loss\_Coupling_$Coupling
	mkdir $dir	
	cp HCN.in submit.sh $dir

	cd $dir
	
	# modify the input parameters
	sed -i "s/in_cavity true 2.01 2 1e-2 100/in_cavity true $Freq 2 $Coupling $Loss/" HCN.in
	sed -i "s/HCN_rtdft/ESC_$dir/" submit.sh
	
	sbatch submit.sh
	cd ..

done    
done
done
