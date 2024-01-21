#!/bin/bash

for dis in 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4
do


dir=d_$dis
echo $dir
mkdir $dir
cd $dir

echo "Au-H dis=$dis ..."

cp ../submit.sh .

RH=$(echo "print($dis+2.222000)" | python3)
RAu1=$(echo "print(2.222000+$dis*2)" | python3)
RAu2=$(echo "print(2.222000+$dis*2+2.222)" | python3)
RAu3=$(echo "print(2.222000+$dis*2+2.222)" | python3)
RAu4=$(echo "print(2.222000+$dis*2+4.633)" | python3)
RAu5=$(echo "print(2.222000+$dis*2+4.537)" | python3)
RAu6=$(echo "print(2.222000+$dis*2+4.542)" | python3)

read -r -d '' str << EOM
\$molecule
+1 1
H      0.000000    0.000000    $RH
Au     0.000000    0.000000    2.222000
Au    -1.392000    0.000000    0.000000
Au     1.393000    0.000000    0.000000
Au     0.000000    0.000000   -2.411000
Au     2.626000    0.000000   -2.315000
Au    -2.614000    0.000000   -2.320000
Au     0.000000    0.000000    $RAu1
Au    -1.392000    0.000000    $RAu2
Au     1.393000    0.000000    $RAu3
Au     0.000000    0.000000    $RAu4
Au     2.626000    0.000000    $RAu5
Au    -2.614000    0.000000    $RAu6
EOM

# add the last line
str="$str\n\$end"

printf "$str" | cat - ../template.in > au12h.in

# submit jobs
sbatch submit.sh

cd ..

done

