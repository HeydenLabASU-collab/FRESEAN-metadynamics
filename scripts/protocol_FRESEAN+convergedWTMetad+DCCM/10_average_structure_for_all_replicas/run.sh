#!/bin/bash
#SBATCH -p htc
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:30                  # wall time (D-HH:MM)
#SBATCH -J REWEIGHT

# Load python enviornment
module load mamba/latest
source activate KEMP_ML

# Input variables
nReplicas=20
path=.. 
inpRefPDB=$path/03-CG/ref.pdb
individual_average_directory=../09_average_structure_for_each_replica
nAtoms=$(( $(wc -l ${inpRefPDB} | awk '{print $1}') - 2 ))
inpWeightedAverageStruct=average-test.dat
outWeightedAverageStruct=average_structure.dat

# Compile f90 code with gfortran
# NOTE: Move .f90 code to separate dir? Compile in make?
gfortran AVERAGE-weighted-pdb.f90 -o AVERAGE-weighted-pdb.out

cat << STOP >& average-reweight-pdb.inp
${outWeightedAverageStruct}
1 ${nAtoms}
${nReplicas}
${individual_average_directory}/template_pdb_0/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_1/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_2/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_3/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_4/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_5/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_6/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_7/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_8/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_9/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_10/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_11/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_12/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_13/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_14/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_15/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_16/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_17/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_18/${inpWeightedAverageStruct}
${individual_average_directory}/template_pdb_19/${inpWeightedAverageStruct}
STOP

# Execute compiled code with provided input file
./AVERAGE-weighted-pdb.out <average-reweight-pdb.inp >& ensemble.out

python3 f95_to_pdb.py  -ref $inpRefPDB -x average_structure.dat -o ensemble_structure.pdb >& f95_to_pdb.out

source deactivate
