#!/bin/bash

nReplicas=20
path=..
reweightFilename=REWEIGHT-distances

for ((j=0;j<$nReplicas;j++))
do

cp -r template_dccm template_dccm_$j
cd template_dccm_$j
sbatch --job-name=dccm$j.run --output=dccm$j.out --export=path=$path,replica=$j,rwFname=$reweightFilename run.sh
cd ..
done
