#!/bin/bash

for dis in 1.3 1.4 1.5 1.6 1.62 1.64 1.66 1.665 1.670 1.675 1.68 1.685 1.690 1.695 1.7 1.72 1.74 1.76 1.78 1.8 1.9 2.0 2.1 2.2 2.3 2.4
do


dir=d_$dis

echo "$dir"
cat $dir/au12h.out | grep -i "E(NEO-SCF)"

done

