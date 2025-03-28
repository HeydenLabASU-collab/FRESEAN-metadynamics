#!/bin/bash

#SBATCH -p general
#SBATCH -G a100:1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH -t 0-04:00                  # wall time (D-HH:MM)
#SBATCH -J EMEQUI

#BEGIN INPUT
inpTOP=../00-prep/topol.top
inpGRO=../00-prep/prep.gro
next=1
nextDir=../02-MD
gmx=gmx_plumed
#END INPUT

files=(
em.mdp
equi.mdp
${inpTOP}
${inpGRO}
)
for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

#Run energy minimization of the freshly solvated protein
if [ ! -d em ]; then
mkdir em
fi
cd em
$gmx grompp -f ../em.mdp -c ../${inpGRO} -p ../${inpTOP} -r ../${inpGRO} -o em.tpr >& grompp.out
$gmx mdrun -v -deffnm em -cpi -nt 1 -pin on >& mdrun.out
cd ..

#Run equilibration MD simulation with position restraints on heavy protein atoms
if [ ! -d equi ]; then
mkdir equi
fi
cd equi
$gmx grompp -f ../equi.mdp -c ../em/em.gro -p ../${inpTOP} -r ../${inpGRO} -o equi.tpr >& grompp.out
$gmx mdrun -v -deffnm equi -cpi >& mdrun.out
cd ..

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f equi/equi.gro ]; then
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
