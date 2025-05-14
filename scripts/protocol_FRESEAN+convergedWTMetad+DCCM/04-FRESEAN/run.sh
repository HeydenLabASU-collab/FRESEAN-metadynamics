#!/bin/bash

#SBATCH -p htc
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -t 0-04:00                  # wall time (D-HH:MM)
#SBATCH -J FRESEAN

#FFTW is required

#BEGIN INPUT

#input topology in custom format
inpMTOP=../03-CG/topol_prot-cg.mtop

#input GMX trajectory file
inpTRR=../03-CG/sample-NPT_prot_pbc-CG.trr

#input reference file rot trans+rot fitting
inpRefGROcg=../03-CG/ref-cg.gro

#number of frames in trajectory
nFrames=1000000

#length of correlation functions in frames
nCorr=100

#stem name of output files
outName=cg

#trajectory time step (picoseconds)
TRJtimestep=0.020
next=1
nextDir=../05-ModeProj
#END INPUT

files=(
${inpMTOP}
${inpTRR}
${inpRefGROcg}
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

#Generate input parameter file for FRESEAN mode analysis
cat << STOP >& matrix.inp
#fnTop (custom .mtop file)
${inpMTOP}
#fnCrd (CP2K position and velocity files)
${inpTRR}
#fnJob (Job file defining atom groups)
static.job
#nRead (Number of frames in the file)
${nFrames}
#analysisInterval (Analysis should be performed every X frames)
1
#fnRef (Reference file generated from genREF.exe)
${inpRefGROcg}
#alignGrp
0
#analyzeGrp
0
#wrap
0
#nCorr (correlation time in trajecotry steps)
${nCorr}
#winSigma (wavenumbers cm^-1)
10.0
#binaryMatrix
1
#doGenModes
0
#convergence (only for generalized normal modes)
1.0e-5
#maxIter (only for generalized normal modes)
100
#fnOut (Output file appender)
${outName}
STOP
#Perform FRESEAN mode analysis for input trajectory
#Here: we use 48 CPU cores --> see SBATCH flag: -c
export OMP_NUM_THREADS=48
fresean matrix -f matrix.inp >& matrix.out

#Diagonalize velocity cross correlation matrix at all sampled frequencies
fresean eigen -m matrix_${outName}.mmat -n $nCorr >& eigen.out

#Generate input parameter file to extract FRESEAN modes sampled at zero frequency
cat << STOP >& extract.inp
#fnEigVec (file containing eigenvectors in binary format)
evec_matrix_${outName}.mmat
#extractMode ( Mode 0 -> freqSel is in wavenumbers; Mode 1 -> freqSel is matrix index )
1
#freqSel (extract mode 0 -> Integer frequency in cm^-1;  mode 1 -> frequency index)
1
#trrFreq
${TRJtimestep}
#modeStart (First vibrational mode to read)
1
#modeEnd (Last vibrational mode to read)
30
#fnOut ( Output file appender )
cg.xyz
STOP
#Extract zero frequency modes
fresean extract -f extract.inp >& extract.out

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f evec_freq1_mode1-30_cg.xyz ]; then
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
