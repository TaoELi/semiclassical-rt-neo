#!/bin/bash

for dt in 1.0 # 0.40 0.80 2.0 8.0 16.0
do
    dir=dt_$dt
    mkdir $dir
    cp HCN.in submit.sh $dir
    cd $dir
    echo "Entering $dir"
    # calculate maxiter according to dt [tmax = 100 fs = 4134.137 a.u.]
    maxiter=$(echo "print(int(4134.137/$dt*10))" | python)
    echo "MaxIter for dt = $dt is $maxiter"
    # calculate field amplitude according to dt (Amp * dt  = 4e-4)
    amp=$(echo "print(4e-4/$dt)" | python)
    echo "Amp for delta pulse is $amp"

    sed -i "s/dt 0.04/dt $dt/" HCN.in
    sed -i "s/field_amp 1e-2/field_amp $amp/" HCN.in
    sed -i "s/maxiter 103353/maxiter $maxiter/" HCN.in
    sed -i "s/BO_pulse_0.04/BO_pulse_$dt/" submit.sh
    
    # submit job
    sbatch submit.sh
    cd .. 
done
