#!/bin/bash

for i in {0..20}
do
    input="oHBA_$i.txt"
    echo $input
    input_new="$i.tmp"
    paste elements.txt $input > $input_new
    
    # create a Qchem input file
    cp oHBA_rgs.in.template oHBA_rgs_$i.in
    sed -i "2r $input_new" oHBA_rgs_$i.in

    rm $input_new
done
