#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-01:00                  # wall time (D-HH:MM)
#SBATCH -J PROJ

# FFTW is required
# Python3 with numpy is required

module load mamba/latest
source activate FRESEAN-python-env

############################
# CHECK NUMPY INSTALLATION #
############################

contains_numpy=`pip list | awk 'BEGIN {i=0} $1 ~ "numpy" {i=1} END {print i}'`

if [ ${contains_numpy} != 1 ];then
echo "You are missing numpy."
exit
fi


############################
#        BEGIN INPUT       #
############################

#all-atom reference structure (must match CG reference)
inpRefPDBaa=../03-CG/ref.pdb

#coarse-grained zero-frequency eigenvectors
inpModesXYZ=../04-FRESEAN/evec_freq1_mode1-30_cg.xyz

# 20 ns unbiased trajectory
inpTRR=../02-MD/sample-NPT_prot_pbc.trr

# If next=1 -> start next step automatically
# If next=0 -> terminate after this step
next=1
nextDir=../06-ModeProj

# If ERROR_FLAG=0 -> no error
ERROR_FLAG=0

############################
#   CHECK FILE EXISTANCE   #
############################

files=(
${inpRefPDBaa}
${inpModesXYZ}
backmap.py
)

for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

###########################
#  BEGIN MODE CONVERSION  #
###########################

cp ${inpRefPDBaa} ref.pdb
cp ${inpModesXYZ} input-cg-modes.xyz

reflinecount=`cat ref.pdb | awk 'BEGIN {natoms=0} $1 ~ "ATOM" || $1 ~ "HETATM" {natoms++} END {print natoms} '`
echo "Found ${reflinecount} atoms in ref.pdb"

# Convert FRESEAN modes 7 & 8 at zero frequency
minmode=1
maxmode=30

for ((mode=${minmode}; mode<=${maxmode}; mode++))
do
python backmap.py $mode >& prep_plumed_mode${mode}.out
modelinecount=`cat evec_${mode}_aa_backmapped.xyz | awk 'BEGIN {natoms=0} $1 ~ "X" {natoms++} END {print natoms}'`
echo "Found ${modelinecount} atoms in evec_${mode}_aa_backmapped.xyz"
if [ ${modelinecount} != ${reflinecount} ]; then
        echo "Mode ${mode} (evec_${mode}_aa_backmapped.xyz) is incomplete."
        ERROR_FLAG=1
fi

done

if [ ${ERROR_FLAG} == 1 ];then
	exit
fi

for i in $(seq $minmode $maxmode); do
    cat evec_${i}_aa_backmapped.xyz
done > evec_aa_${minmode}-${maxmode}_backmapped.xyz

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
    if [ -f evec_aa_${minmode}-${maxmode}_backmapped.xyz ]; then
        if [ -d ${nextDir} ]; then
            curDir=`pwd`
            cd ${nextDir}
            if [ -f run.sh ]; then
                sbatch run.sh
            fi
            cd ${curDir}
        fi
    fi
fi
