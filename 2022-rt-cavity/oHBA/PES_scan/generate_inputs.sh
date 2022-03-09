#!/bin/bash

for ratio in $(seq 0 0.01 0.5)
do
    filename=oHBA_rgs_$ratio
    
    # create qchem input file
    ./create_qchem_in.sh $ratio $filename.in
    
    # create submit file
    cp submit.sh submit_$ratio.sh
    sed -i "s/qchem oHBA_rgs.in oHBA_rgs.out/qchem $filename.in $filename.out/" submit_$ratio.sh
    
    # submit job
    sbatch submit_$ratio.sh

done
