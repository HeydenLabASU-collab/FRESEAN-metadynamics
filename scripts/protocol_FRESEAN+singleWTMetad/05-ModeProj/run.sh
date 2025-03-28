#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-01:00                  # wall time (D-HH:MM)
#SBATCH -J PROJ

# FFTW is required
# Python3 with numpy is required

#BEGIN INPUT
#all-atom reference structure (must match CG reference)
inpRefPDBaa=../03-CG/ref.pdb
inpModesXYZ=../04-FRESEAN/evec_freq1_mode1-30_cg.xyz
inpTRR=../02-MD/sample-NPT_prot_pbc.trr
next=1
nextDir=../06-metadyn
#END INPUT

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

#Convert FRESEAN modes 7 & 8 at zero frequency
cp ${inpRefPDBaa} ref.pdb
cp ${inpModesXYZ} input-cg-modes.xyz
mode=7
python prep_plumed.py $mode >& prep_plumed_mode${mode}.out
mode=8
python prep_plumed.py $mode >& prep_plumed_mode${mode}.out


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
