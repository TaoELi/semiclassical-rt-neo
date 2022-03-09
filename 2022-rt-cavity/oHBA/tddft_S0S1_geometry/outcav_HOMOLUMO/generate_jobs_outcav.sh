#!/bin/bash

for CavFreq in 3.2959
do
for CavLossTime in 1e8
do
for Coupling in 0e-4
do
for Amp in 0e-2
do
for i in {0..20}
do

	#dir=Geom_$i\_CavFreq_$CavFreq\_Coupling_$Coupling\_Loss_$CavLossTime\_Amp_$Amp
	dir=Geom_$i
	echo "generating $dir"
	mkdir $dir	
	cp oHBA_rgs_$i.in submit.sh $dir

	cd $dir
	
	# modify the input parameters
	sed -i "s/field_amp 0.000/field_amp $Amp/" oHBA_rgs_$i.in
	sed -i "s/HCN_rtdft/s0s1_$dir/" submit.sh
	#sed -i "s/in_cavity false 2.01 2 1e-2 100/in_cavity true $CavFreq 1 $Coupling $CavLossTime/" oHBA_rgs_$i.in
	
	sed -i "s/oHBA_rgs.in/oHBA_rgs_$i.in/" submit.sh
	sbatch submit.sh
	cd ..

done
done
done
done
done
