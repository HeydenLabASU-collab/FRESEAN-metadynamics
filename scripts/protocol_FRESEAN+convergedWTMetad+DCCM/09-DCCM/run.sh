#!/bin/bash
#SBATCH -p htc
#SBATCH -q public
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-04:00
#SBATCH -J DCCM

module load mamba/latest
source deactivate
source activate mda

cd ctypes
./compile.sh
cd ..

python dccm.py