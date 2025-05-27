#!/bin/bash

for ((j=0;j<20;j++))
do
cp -r single_rw reweight_$j
sed -i "7s/0/$j/" reweight_$j/startme.sh
sed -i "14s/0/$j/" reweight_$j/startme.sh
cd reweight_$j
sbatch startme.sh
cd ..
done
