#!/bin/bash

for dir in d_*/
do

echo "$dir"
cd $dir

cat au12h.out.esp | tail -n14001 > esp.txt

cd ..
done

