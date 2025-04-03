#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-08:00                  # wall time (D-HH:MM)
#SBATCH -J COARSE

# FFTW is required

#BEGIN INPUT
gmx=gmx_plumed

#GMX topology for protein
inpTOPprot=../00-prep/topol_prot.top

#GMX files from MD sampling
inpTPR=../02-MD/sample-NPT.tpr
inpTRR=../02-MD/sample-NPT.trr

#output trajectory for protein with PBC fixes
outTRRprotAA=../02-MD/sample-NPT_prot_pbc.trr

#output trajectory for protein in CG representation
outTRRprotCG=sample-NPT_prot_pbc-CG.trr

#GMX group number for protein (gmx trjconv)
outGrp=1

#time between frames in input trajectory
TRJtimestep=0.020

#number of steps in input trajectory
TRJframes=1000000
next=1
nextDir=../04-FRESEAN
#END INPUT

files=(
${inpTOPprot}
${inpTPR}
${inpTRR}
static.job
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Generate a all-atom topology file (custom format) for the protein
fresean mtop << STOP >& mtop.out
${inpTOPprot}
${inpTOPprot}
topol_prot-aa.mtop
STOP

#Make the protein whole (fix PBC jumps) and write out protein-only trajectory
$gmx trjconv -s ${inpTPR} -f ${inpTRR} -o ${outTRRprotAA} -pbc mol << STOP >& trjconv1.out
${outGrp}
STOP

#Generate all-atom reference structure for rot+trans fitting
$gmx trjconv -s ${inpTPR} -f ${outTRRprotAA} -o ref.gro -e 0.0 << STOP >& trjconv2.out
${outGrp}
STOP
#some postprocessing for PLUMED compatibility of ref.pdb
$gmx editconf -f ref.gro -o ref.pdb
tail -n +5 ref.pdb > tmp.pdb
head -n -2 tmp.pdb > ref.pdb
sed -i '1i REMARK TYPE=OPTIMAL' ref.pdb
sed -i '$aEND' ref.pdb
rm tmp.pdb

#Generate input parameter file for coarse-graining
#WARNING: update expected soon
cat << STOP >& coarse.inp
#fnTop
topol_prot-aa.mtop
#fnCrd
${outTRRprotAA}
#fnVel (if format xyz,crd,dcd)
#fnJob
static.job
#grp
0
#nRead
${TRJframes}
#analysisInterval
1
#fnOutTraj
tmptraj.gro
#fnOutTopol
topol_prot-cg.mtop
STOP
#Apply coarse-graining to all-atom protein trajectory
fresean coarse -f coarse.inp >& coarse.out
#Convert coarse-grained trajectory from gro (ASCII text) to trr (binary) format
$gmx trjconv -s tmptraj.gro -f tmptraj.gro -timestep ${TRJtimestep} -o ${outTRRprotCG} << STOP >& trjconv3.out
0
STOP

#Generate coarse-grained reference structure for rot+trans fitting
nAtoms=`head -n 2 tmptraj.gro | tail -n 1`
nLines=${nAtoms}
((nLines +=3))
head -n ${nLines} tmptraj.gro >& ref-cg.gro
rm tmptraj.gro

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f ${outTRRprotCG} ]; then
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
