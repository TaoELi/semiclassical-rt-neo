#!/bin/bash

for N in 16 20 24 26 28 30 32 34
do

for dis in 0.0
do


dir=N_$N\_d1_$dis\_excite_all_z
echo $dir
mkdir $dir
cd $dir

echo "N=$N Na atoms, Na-HCN dis=$dis ..."

cp ../submit.sh .

RH=$(echo "print($dis)" | python3)
RF=$(echo "print($dis-0.9245837)" | python3)

read -r -d '' str << EOM
\$molecule
0 1
H        $RH    0.0000000000    0.0000000000
F        $RF    0.0000000000    0.0000000000 
EOM

# create the geometry for Na
for (( idx=1; idx<=$N; idx++ ))
do
#echo "$N, $idx"
disNa=$(echo "print(($idx-1)*2.55 + 2.2)" | python3)
str="$str \nAu $disNa 0.0 0.0"
#echo "$str"
done

# add the last line
str="$str\n\$end"

printf "$str" | cat - ../template_excite_all.in > hfau2n.in

sed -i "s/field_type delta/field_type delta\nfield_direction 2/" hfau2n.in

# change the number of cpus if the number is very small
if [ 8 -gt $N ] 
then
sed -i "s/32/8/g" submit.sh
fi

# submit jobs
sed -i "s/job-name=au/job-name=au_excite_all_z/" submit.sh
sbatch submit.sh

cd ..

done

done

