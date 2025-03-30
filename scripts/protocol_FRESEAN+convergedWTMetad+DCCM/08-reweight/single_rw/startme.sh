#!/bin/bash
#SBATCH -p lightwork
#SBATCH -q public
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:30
#SBATCH -J REWEIGHT_0

#BEGIN INPUT

gmx=gmx_plumed
#GMX topology for protein
workingDir=../..
i=${replica}
inpTPRprot=$workingDir/07-metadyn/metadyn_$i/metadyn_prot.tpr
inpTRRprot=$workingDir/07-metadyn/metadyn_$i/metadyn_prot_pbc.trr
inpPlumedPDB=$workingDir/05-ModeProj/plumed-mode-input.pdb
inpPlumedHills=$workingDir/07-metadyn/metadyn_$i/plumed-mode-metadyn.hills

echo "Input TPR $inpTPRprot"
echo "Input TRR $inpTRRprot"
#END INPUT

files=(
${inpTPRprot}
${inpTRRprot}
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Copy PLUMED input file with reference structure and FRESEAN modes to standardized file name
cp ${inpPlumedPDB} plumed-mode-input.pdb
#Copy PLUMED hills file from metadynamics run as input for reweighting
cp ${inpPlumedHills} plumed-mode-metadyn.hills

$gmx mdrun -s ${inpTPRprot} -nsteps 1 -plumed ../plumed-mass+charge.dat >& plumed-mc.out

# kT in kJ/mol can be determined by running `plumed kT --temp 300`
kt=2.494339

#Generate groups.ndx file from metadynamics gro output
echo -e "a CA\nname 13 C-alpha\nq" | gmx make_ndx -f confout.gro -o groups.ndx >& make-ndx.out

#Generate reweighted histogram as a function of new CVs
plumed driver --mf_trr ${inpTRRprot} --plumed ../plumed-reweight-CV.dat --kt $kt --mc mass+charge.dat > reweight.out
