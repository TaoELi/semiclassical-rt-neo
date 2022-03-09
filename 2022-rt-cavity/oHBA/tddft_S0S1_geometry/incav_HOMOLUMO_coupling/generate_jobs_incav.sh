#!/bin/bash

CavFreqArray=(3.295 3.329 3.363 3.396 3.429 3.460 3.491 3.522 3.553 3.582 3.611 3.640 3.668 3.696 3.725 3.753 3.780 3.808 3.836 3.863 3.890)

for CavLossTime in 1e8
do
for Coupling in 6e-3 #1e-3 2e-3 3e-3 4e-3 5e-3 6e-3
do
for Amp in 0e-2
do
for i in 0 20
do

	CavFreq=${CavFreqArray[$i]}
	echo "CavFreq=$CavFreq"
	#dir=Geom_$i\_CavFreq_$CavFreq\_Coupling_$Coupling\_Loss_$CavLossTime\_Amp_$Amp
	dir=Geom_$i\_Coupling_$Coupling
	echo "generating $dir"
	mkdir $dir	
	cp oHBA_rgs_$i.in submit.sh $dir

	cd $dir
	
	# modify the input parameters
	sed -i "s/field_amp 0.000/field_amp $Amp/" oHBA_rgs_$i.in
	sed -i "s/HCN_rtdft/incav_$dir/" submit.sh
	sed -i "s/in_cavity false 2.01 2 1e-2 100/in_cavity true $CavFreq 1 $Coupling $CavLossTime/" oHBA_rgs_$i.in
	
	sed -i "s/oHBA_rgs.in/oHBA_rgs_$i.in/" submit.sh
	sbatch submit.sh
	cd ..

done
done
done
done
