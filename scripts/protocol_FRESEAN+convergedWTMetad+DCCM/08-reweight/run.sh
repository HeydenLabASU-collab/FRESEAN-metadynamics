#!/bin/bash

nReplicas=20
for (( i=0; i<${nReplicas}; i++ ))
do
cp -r single_rw reweight_$i
cd reweight_$i
echo "Run replica $i"
sbatch --job-name=REWEIGHT_REPLICA_${i}.run --export=replica=${i} startme.sh
cd ..
done
