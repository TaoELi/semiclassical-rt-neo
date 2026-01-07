#!/bin/bash

dir=$1

echo "working in dir $dir"

cd $dir

echo "collecting dipole info"
sed 's/^[ \t]*//' *.out | awk '/full coordinates/{getline; print}' > x.tmp
wait
sed 's/^[ \t]*//' *.out | awk '/full coordinates/{getline; getline; print}' > y.tmp
wait
sed 's/^[ \t]*//' *.out | awk '/full coordinates/{getline; getline; getline; print}' > z.tmp
wait

cat *.out | grep -i "mu_n is" | awk '{print $3 " " $6 " " $7 " " $8}' > xn.txt
cat *.out | grep -i "e_conserved is" | awk '{print $3 " " $6}' > e_tot.txt
cat *.out | grep -i "e_kin(quantum) is" | awk '{print $3 " " $6}' > e_kq.txt
cat *.out | grep -i "e_kin(qmclcross) is" | awk '{print $3 " " $6}' > e_cr.txt
cat *.out | grep -i "e_kin(qmqm) is" | awk '{print $3 " " $6}' > e_qm.txt
cat *.out | grep -i "e_constraint_pbc is" | awk '{print $3 " " $6}' > e_pbc.txt
cat *.out | grep -B2 -I "Time step" |  grep -i "Tau rescaling" | awk '{print $12}' > tau_dynamics.txt
grep -B4 -i "Time step " *.out| grep -i "ftot =" | awk '{print $14}' > ftot_dynamics.txt
grep -B4 -i "Time step " *.out| grep -i "ftot - ntot" | awk '{print $6}' > ftot2_dynamics.txt

n=$(wc -l x.tmp | awk '{print $1}')
echo "the total length of the simulation is $n"

#dt=4
#tmax=$(python -c "print($dt*($n-1))")
#echo "tmax is $tmax"

#seq 0.0 $dt $tmax > t_series.tmp
#wait

grep -i "full coordinates" *.out | awk '{print $4}' > t_series.tmp

paste t_series.tmp x.tmp > x.txt
wait
paste t_series.tmp y.tmp > y.txt
wait
paste t_series.tmp z.tmp > z.txt
wait

rm t_series.tmp x.tmp y.tmp z.tmp

cd ..
