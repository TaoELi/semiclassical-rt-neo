#!/bin/bash

for CavFreq in 3.2959
do
for CavLossTime in 1e8
do
for Coupling in 0e-4
do
for Amp in 1e-2
do
for i in {0..20}
do

	dir=Geom_$i\_LR
	echo "generating $dir"
	mkdir $dir	
	cp oHBA_rgs_$i.in submit.sh $dir

	cd $dir
	
	# modify the input parameters
	sed -i "s/field_amp 0.000/field_amp $Amp/" oHBA_rgs_$i.in
	sed -i "s/electronic_HOMO2LUMO true/electronic_HOMO2LUMO false/" oHBA_rgs_$i.in
	sed -i "s/HCN_rtdft/LR_$dir/" submit.sh
	
	sed -i "s/oHBA_rgs.in/oHBA_rgs_$i.in/" submit.sh
	sbatch submit.sh
	cd ..

done
done
done
done
done
