#!/bin/bash                                                                     
#SBATCH -p general
#SBATCH -G a100:1
#SBATCH -c 12
#SBATCH -N 1
#SBATCH -t 7-00:00
#SBATCH -J METAD_0

# FFTW IS REQUIRED (PLUMED)

# Working directory containing starting topology
workingDir=../..
inpTOP=$workingDir/00-prep/topol.top
gmx=gmx_plumed

# Starting struct
inpGRO=$workingDir/06-resample/snapshots/state_${replica}.gro

# MetaD CV file
inpPlumedPDB=$workingDir/05-ModeProj/plumed-mode-input.pdb
#GMX group number for protein (gmx trjconv)
outGrp=1
#END INPUT

files=(
${inpTOP}
${inpGRO}
${inpPlumedPDB}
metadyn.mdp
plumed-mode-metadyn.dat
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

if [ ! -f metadyn.tpr ]; then
#Create tpr input file for WT-metadynamics simulation
$gmx grompp -f metadyn.mdp -c ${inpGRO} -p ${inpTOP} -o metadyn.tpr >& grompp.out
fi

if [ ! -f metadyn_prot.tpr ]; then
#Create tpr file with only protein atoms for analysis
$gmx convert-tpr -s metadyn.tpr -o metadyn_prot.tpr << STOP >& tpr-convert.out
${outGrp}
STOP
fi

#Run the WT-metadynamics simulation with GROMACS & PLUMED
$gmx mdrun -v -deffnm metadyn -cpi -plumed plumed-mode-metadyn.dat >& mdrun.out

# kT in kJ/mol can be determined by running `plumed kT --temp 300`
kt=2.494339
#Generate free energy surface in FRESEAN mode space
plumed sum_hills --min -0.02,-0.02 --max 0.02,0.02 --bin 200,200 --hills plumed-mode-metadyn.hills --outfile plumed-mode-metadyn.fes --mintozero --kt 2.494339 --stride 5000

#Write out metadynamics trajectory for only protein atoms after PBC fix
$gmx trjconv -s metadyn.tpr -f metadyn.trr -o metadyn_prot_pbc.trr -pbc mol << STOP >& trjconv.out
${outGrp}
STOP

rm metadyn.trr
