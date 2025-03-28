#!/bin/bash

#SBATCH -p general
#SBATCH -G a100:1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH -t 0-12:00                  # wall time (D-HH:MM)
#SBATCH -J SAMPLE

#BEGIN INPUT
inpTOP=../00-prep/topol.top
inpGRO=../01-em+equi/equi/equi.gro
next=1
nextDir=../03-CG
gmx=gmx_plumed
#END INPUT

files=(
sample-NPT.mdp
${inpTOP}
${inptGRO}
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Run the equilibrium MD simulation to sample vibrational modes
#Here: 1GPU and 16 CPU cores are used --> see SBATCH flags: -G -c
$gmx grompp -f sample-NPT.mdp -c ${inpGRO} -p ${inpTOP} -o sample-NPT.tpr >& grompp.out
$gmx mdrun -v -deffnm sample-NPT -cpi >& mdrun.out

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f sample-NPT.gro ]; then
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
