#!/bin/bash


for dir in rtBO_delta_proton2/
#for dir in rt_delta_proton/
do
    cd $dir
    for dt in 0.40 1.0 4.0 16.0
    do
        pathname="dt_$dt"
        cd $pathname
        echo "working under $dir/$pathname"
        
        filename="HCN.out"
        # get protonic dipole moment traj
        grep -i "mu_n is" $filename | awk '{print $6 " " $7 " " $8}' > mu_n.tmp
        # get nstep series
        grep -i "mu_n is" $filename | awk '{print $3}' > nstep.tmp
        # multiple nstep by dt to get time series
        nmax=$(wc -l nstep.tmp | awk '{print $1}')
        tmax=$(echo "print(($nmax-1)*$dt)" | python)
        echo "dt=$dt, nstep=$nmax, create time series with tmax=$tmax"
        seq 0 $dt $tmax > tseries.tmp

        paste nstep.tmp tseries.tmp mu_n.tmp > mu_n.txt
        # delete the last line of mu_n.txt
        sed -i '' -e '$ d' mu_n.txt
        rm nstep.tmp tseries.tmp mu_n.tmp

        cd ../
    done
    cd ../
done

