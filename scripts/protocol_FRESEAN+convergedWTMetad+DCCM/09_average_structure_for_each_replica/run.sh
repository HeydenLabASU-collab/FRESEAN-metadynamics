#!/bin/bash

nReplicas=20
path=..
reweightFilename=REWEIGHT-distances
# For each replica
for ((j=0;j<nReplicas;j++))
do

# Copy the template directory to template_pdb_replicaNumber
cp -r template_pdb template_pdb_$j

# DO NOT USE SED TO MODIFY SBATCH .SH SCRIPT JUST USE BUILT-IN EXPORTS

# Generate the average structure for each replica in parallel
cd template_pdb_$j
sbatch --job-name=avg$j.run --output=avg$j.out --export=path=$path,replica=$j,rwFname=$reweightFilename run_generate_weighted_average_structure.sh
cd ..
done
