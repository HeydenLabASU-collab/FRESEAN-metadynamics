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
nextDir=../06-resample

# If ERROR_FLAG=0 -> no error
ERROR_FLAG=0

############################
#   CHECK FILE EXISTANCE   #
############################

files=(
${inpRefPDBaa}
${inpModesXYZ}
prep_plumed.py
plumed-mode-projection.dat
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

# Convert FRESEAN modes 7 & 8 at zero frequency
mode=7
python prep_plumed.py $mode >& prep_plumed_mode${mode}.out
mode=8
python prep_plumed.py $mode >& prep_plumed_mode${mode}.out

###########################
#  CHECK PREP COMPLETION  #
###########################

reflinecount=`cat ref.pdb | awk 'BEGIN {natoms=0} $1 ~ "ATOM" || $1 ~ "HETATM" {natoms++} END {print natoms} '`
mode7linecount=`cat evec_7_aa_scaled.pdb | awk 'BEGIN {natoms=0} $1 ~ "ATOM" || $1 ~ "HETATM" {natoms++} END {print natoms}'`
mode8linecount=`cat evec_8_aa_scaled.pdb | awk 'BEGIN {natoms=0} $1 ~ "ATOM" || $1 ~ "HETATM" {natoms++} END {print natoms}'`
echo "Found ${reflinecount} atoms in ref.pdb"
echo "Found ${mode7linecount} atoms in evec_7_aa_scaled.pdb"
echo "Found ${mode8linecount} atoms in evec_8_aa_scaled.pdb"

if [ ${mode7linecount} != ${reflinecount} ]; then
	echo "Mode 7 (evec_7_aa_scaled.pdb) is incomplete."
	ERROR_FLAG=1
fi

if [ ${mode8linecount} != ${reflinecount} ]; then
	echo "Mode 8 (evec_8_aa_scaled.pdb) is incomplete."
	ERROR_FLAG=1
fi

if [ ${ERROR_FLAG} == 1 ];then
	exit
fi

#Combine extracted modes with all-atom reference input PLUMED input
cat ref.pdb evec_7_aa_scaled.pdb evec_8_aa_scaled.pdb > plumed-mode-input.pdb
rm ref.pdb input-cg-modes.xyz evec_7_aa_scaled.pdb evec_8_aa_scaled.pdb

#Analyze displacement fluctuations along FRESEAN modes
plumed driver --mf_trr ${inpTRR} --plumed plumed-mode-projection.dat --kt 2.494339 > plumed-driver.out

#number of histogram bins
bins=200

#print standard deviations of displacement fluctuations along FRSEAN modes
python standard-deviation.py plumed-mode-projection.out $bins >& standard-deviation.out

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f plumed-mode-projection.dat ]; then
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
