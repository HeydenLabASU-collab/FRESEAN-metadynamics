#!/bin/bash
#SBATCH -p general
#SBATCH -q public
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-04:00                  # wall time (D-HH:MM)                       
#SBATCH -J REWEIGHT_REPLICA_0

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

inpRefPDB=$workingDir/03-CG/ref.pdb
rw_file=$workingDir/08-reweight/reweight_$replica_id/$reweightFilename

weighted_average_structure=../../10_average_structure_for_all_replicas/ensemble_structure.pdb
inpXTCprot=../../09_average_structure_for_each_replica/template_pdb_$replica_id/metadyn_prot_fitted.xtc

echo "Input weighted_average_structure $weighted_average_structure"
echo "Input XTC $inpXTCprot"
#END INPUT

files=(
${weighted_average_structure}
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

# Fit the trajectory to the weighted average structure
$gmx trjconv -s ${weighted_average_structure} -f ${inpXTCprot} -o metadyn_prot_fit_to_avg.xtc -fit rot+trans << STOP >& trjconv.out
${outGrp}
${outGrp}
STOP

python3 dccm_reweight.py -ref $inpRefPDB -trr metadyn_prot_fit_to_avg.xtc -cv ${rw_file} -av $weighted_average_structure -o dccm.dat >& dccm.out

source deactivate
