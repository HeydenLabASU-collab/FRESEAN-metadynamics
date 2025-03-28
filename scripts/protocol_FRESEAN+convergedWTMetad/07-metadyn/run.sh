#!/bin/bash

# For each replica
for ((i=0;i<20;i++));
do

# Copy template
cp -r single_metad metadyn_$i
cd metadyn_$i

# Change the line that points to what snapshot to use
echo "Run replica $i"
sed -i '7s/\(^\|[^a-zA-Z0-9]\)0\([^a-zA-Z0-9]\|$\)/\1'"$i"'\2/g' startme.sh
sed -i '15s/\(^\|[^a-zA-Z0-9]\)0\([^a-zA-Z0-9]\|$\)/\1'"$i"'\2/g' startme.sh
sbatch startme.sh
cd ..
done
