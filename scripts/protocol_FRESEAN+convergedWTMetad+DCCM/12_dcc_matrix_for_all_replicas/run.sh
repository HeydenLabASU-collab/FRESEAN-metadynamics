#!/bin/bash
#SBATCH -p htc
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:30                  # wall time (D-HH:MM)
#SBATCH -J REWEIGHT

# Input variables
nReplicas=20
path=.. 
inpRefPDB=$path/03-CG/ref.pdb
individual_average_directory=../11_dcc_matrix_for_each_replica
nCA=$(awk '$3 == "CA" { count++ } END { print count }' $inpRefPDB)
inpDCCM=dccm.dat
outDCCM=dccm-reweight.dat

# Compile f90 code with gfortran
# NOTE: Move .f90 code to separate dir? Compile in make?
gfortran AVERAGE-dccm.f90 -o AVERAGE-dccm.out

cat << STOP >& average-dccm.inp
${outDCCM}
1 ${nCA}
${nReplicas}
${individual_average_directory}/template_dccm_0/${inpDCCM}
${individual_average_directory}/template_dccm_1/${inpDCCM}
${individual_average_directory}/template_dccm_2/${inpDCCM}
${individual_average_directory}/template_dccm_3/${inpDCCM}
${individual_average_directory}/template_dccm_4/${inpDCCM}
${individual_average_directory}/template_dccm_5/${inpDCCM}
${individual_average_directory}/template_dccm_6/${inpDCCM}
${individual_average_directory}/template_dccm_7/${inpDCCM}
${individual_average_directory}/template_dccm_8/${inpDCCM}
${individual_average_directory}/template_dccm_9/${inpDCCM}
${individual_average_directory}/template_dccm_10/${inpDCCM}
${individual_average_directory}/template_dccm_11/${inpDCCM}
${individual_average_directory}/template_dccm_12/${inpDCCM}
${individual_average_directory}/template_dccm_13/${inpDCCM}
${individual_average_directory}/template_dccm_14/${inpDCCM}
${individual_average_directory}/template_dccm_15/${inpDCCM}
${individual_average_directory}/template_dccm_16/${inpDCCM}
${individual_average_directory}/template_dccm_17/${inpDCCM}
${individual_average_directory}/template_dccm_18/${inpDCCM}
${individual_average_directory}/template_dccm_19/${inpDCCM}
STOP

# Execute compiled code with provided input file
./AVERAGE-dccm.out <average-dccm.inp >& average-dccm-driver.out
