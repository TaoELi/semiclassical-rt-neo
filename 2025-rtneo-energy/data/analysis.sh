#!/bin/bash

for dir in pgrad_0_mv_0_rescf_0_dtsm/ 
do

./get_full_coordinates.sh $dir
python generate_coordinate_oHBA.py $dir
python generate_coordinate_oHBA_sub.py $dir
python generate_coordinate_oHBA_bc.py $dir

done
