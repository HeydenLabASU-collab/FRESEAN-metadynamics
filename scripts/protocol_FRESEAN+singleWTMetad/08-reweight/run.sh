#!/bin/bash
#SBATCH -p htc
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-04:00
#SBATCH -J REWEIGHT

#BEGIN INPUT

gmx=gmx_plumed
#GMX topology for protein
inpTPRprot=../07-metadyn/metadyn_prot.tpr
inpTRRprot=../07-metadyn/metadyn_prot_pbc.trr
inpPlumedPDB=../07-metadyn/plumed-mode-input.pdb
inpPlumedHills=../07-metadyn/plumed-mode-metadyn.hills
#END INPUT

files=(
${inpTPRprot}
${inpTRRprot}
plumed-mass+charge.dat
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Extract atom masses and atom charges from tpr file with only protein atoms
$gmx mdrun -s ${inpTPRprot} -nsteps 1 -plumed plumed-mass+charge.dat >& plumed-mc.out

#Generate GROMACS index (ndx) file with standard groups
#Modify as needed, e.g., if specific groups are used to define CVs for PLUMED
$gmx make_ndx -f ${inpTPRprot} -o groups.ndx << STOP >& make_ndx.out
q
STOP

#Copy PLUMED input file with reference structure and FRESEAN modes to standardized file name
cp ${inpPlumedPDB} plumed-mode-input.pdb
#Copy PLUMED hills file from metadynamics run as input for reweighting
cp ${inpPlumedHills} plumed-mode-metadyn.hills

# kT in kJ/mol can be determined by running `plumed kT --temp 300`
kt=2.494339
#Generate reweighted histogram as a function of new CVs
plumed driver --mf_trr ${inpTRRprot} --plumed plumed-reweight-CV.dat --kt $kt --mc mass+charge.dat > reweight.out
