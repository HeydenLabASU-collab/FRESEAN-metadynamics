#!/bin/bash
#SBATCH -p htc
#SBATCH -q public
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:30                  # wall time (D-HH:MM)                       
#SBATCH -J DCCM_STEP1_0

# Module Loading
module load mamba/latest
source activate KEMP_ML

# Constants
gmx=gmx_plumed
outGrp=1 #protein

# Input Variables from SBATCH command
replica_id=$replica
workingDir=$path
reweightFilename=$rwFname

# Input Variables
inpTPRprot=$workingDir/07-metadyn/metadyn_0/metadyn_prot.tpr
inpTRRprot=$workingDir/07-metadyn/metadyn_$replica_id/metadyn_prot_pbc.trr
inpRefPDB=$workingDir/03-CG/ref.pdb 
rw_file=$workingDir/08-reweight/reweight_$replica_id/$reweightFilename

echo "Input TPR $inpTPRprot"
echo "Input TRR $inpTRRprot"
echo "Input Reference PDB $inpRefPDB"
echo "Input Reweight Filepath $rw_file"

#END INPUT

files=(
${inpTPRprot}
${inpTRRprot}
${inpRefPDB}
${rw_file}
)

for file in ${files[@]}
do
if [ ! -f ${file} ]; then
echo "-could not find file ${file} in current directory"
echo "-exiting"
exit
fi
done

# Generate protein trajectory fit to reference pdb
$gmx trjconv -s ${inpRefPDB} -f ${inpTRRprot} -o metadyn_prot_fitted.xtc -fit rot+trans << STOP >& trjconv.out
${outGrp}
${outGrp}
STOP

# Generate the weighted average structure
python3 generate_weighted_average_structure.py -ref $inpRefPDB -trr metadyn_prot_fitted.xtc -cv ${rw_file} -o average-test.dat >& average.out

# Deactivate python enviornment
source deactivate
