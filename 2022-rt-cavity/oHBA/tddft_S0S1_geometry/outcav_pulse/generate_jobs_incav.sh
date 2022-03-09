#!/bin/bash

CavFreqArray=(3.295 3.329 3.363 3.396 3.429 3.460 3.491 3.522 3.553 3.582 3.611 3.640 3.668 3.696 3.725 3.753 3.780 3.808 3.836 3.863 3.890)
LPFreqArray=(1.0    1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0   3.355   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0)

for CavLossTime in 1e8
do
for Coupling in 4e-3
do
for Amp in 2e-2 4e-2 8e-2
do
for i in 10
do

	CavFreq=${CavFreqArray[$i]}
	echo "CavFreq=$CavFreq"
	LP_freq=${LPFreqArray[$i]}
	echo "LP_freq=$LP_freq"
	#dir=Geom_$i\_CavFreq_$CavFreq\_Coupling_$Coupling\_Loss_$CavLossTime\_Amp_$Amp
	dir=Geom_$i\_Amp_$Amp
	echo "generating $dir"
	mkdir $dir	
	cp oHBA_rgs_$i.in submit.sh $dir

	cd $dir
	
	# modify the input parameters
	sed -i "s/field_amp 0.000/field_amp $Amp/" oHBA_rgs_$i.in
	sed -i "s/HCN_rtdft/outcav_$dir/" submit.sh
	#sed -i "s/in_cavity false 2.01 2 1e-2 100/in_cavity true $CavFreq 1 $Coupling $CavLossTime/" oHBA_rgs_$i.in
	
	sed -i "s/oHBA_rgs.in/oHBA_rgs_$i.in/" submit.sh
	
	# we now change the HOMO 2 LUMO transition --> pulse excitation of molecular frequency (at resonance with the cavity frequency)
	sed -i "s/field_type delta/field_type gaussian 0.0 400.0 $CavFreq/" oHBA_rgs_$i.in
 	sed -i "s/field_amp 0.000/field_amp $Amp/" oHBA_rgs_$i.in
	sed -i 's/electronic_HOMO2LUMO true/FIELD_DIRECTION 1/' oHBA_rgs_$i.in

	sbatch submit.sh
	cd ..

done
done
done
done
