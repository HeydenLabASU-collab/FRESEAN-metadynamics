#!/bin/bash
#SBATCH -p general
#SBATCH -G a100:1
#SBATCH -c 12
#SBATCH -N 1
#SBATCH -t 4-00:00
#SBATCH -J METAD

#BEGIN INPUT
inpTOP=../00-prep/topol.top
gmx=gmx_plumed
inpGRO=../01-em+equi/equi/equi.gro
inpPlumedPDB=../06-ModeProj/plumed-mode-input.pdb
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


vals_file="../06-ModeProj/std_dev.txt"
in_file="plumed-mode-metadyn.dat"

# Read first two lines
line1="$(sed -n '1p' "$vals_file")"
line2="$(sed -n '2p' "$vals_file")"

# Escape for sed replacement (handles \, &, and delimiter |)
escape_sed_repl() {
  printf '%s' "$1" | sed 's/[\/&|\\]/\\&/g'
}

r1="$(escape_sed_repl "$line1")"
r2="$(escape_sed_repl "$line2")"

# Replace XXX -> line1 and YYY -> line2
sed -i -e "s|XXX|$r1|g" -e "s|YYY|$r2|g" "$in_file"

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
plumed sum_hills --hills plumed-mode-metadyn.hills --outfile plumed-mode-metadyn.fes  --mintozero --kt $kt --stride 5000

#Write out metadynamics trajectory for only protein atoms after PBC fix
$gmx trjconv -s metadyn.tpr -f metadyn.trr -o metadyn_prot_pbc.trr -pbc mol << STOP >& trjconv.out
${outGrp}
STOP
