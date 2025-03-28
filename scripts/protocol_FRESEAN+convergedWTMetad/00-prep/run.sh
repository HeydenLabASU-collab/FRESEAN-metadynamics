#!/bin/bash

#BEGIN INPUT
gmx=gmx_plumed
pdb=1hel.pdb
ff=6 #Option 6 corresponds to AMBER99sb-ILDN
wat=1 #Option 1 corresponds to tip3p
box=9.3 #Default box size
salt=0.15 #Default salt concentration
ionGrp=13
#END INPUT

#Generate first topology file (protein + crystal water/ions)
$gmx pdb2gmx -f ${pdb} -p topol_prot.top -o prot.gro << STOP
${ff}
${wat}
STOP

#Define simulation box for preiodic boundary conditions (PBC)
$gmx editconf -f prot.gro -box ${box} -o box.gro

#Make copy of protein topology file and delete crystal water/ions from original
cp topol_prot.top topol.top
awk '{if($1!="SOL"&&$1!="NA"&&$1!="CL") {printf("%s\n",$0);}}' topol.top >& topol_prot.top

#Solvate the protein and update new topology file (keep crystal water/ions)
$gmx solvate -cp box.gro -cs -p -o solv.gro
rm \#topol.top.1\#

#Add ions to the solution (requires *.tpr input file)
$gmx grompp -f em.mdp -c solv.gro -p -o tmp.tpr -maxwarn 1
$gmx genion -s tmp.tpr -p -o prep.gro -neutral -conc ${salt} << STOP
${ionGrp}
STOP
rm \#topol.top.1\#
