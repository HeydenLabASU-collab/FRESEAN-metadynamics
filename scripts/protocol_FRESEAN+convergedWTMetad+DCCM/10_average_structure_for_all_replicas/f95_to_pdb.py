import MDAnalysis as mda
import numpy as np
import math
import argparse

def read_coordinate_data(coordinate_file):
    coordinates = np.loadtxt(coordinate_file) 
    return coordinates

#Parse arguments for reference pdb and the trajectory
parser = argparse.ArgumentParser()
parser.add_argument("-ref", "--reference", type=str, required=True)
parser.add_argument("-x", "--coordinates", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

reference_pdb_file = args.reference
coordinate_file = args.coordinates
output_pdb_file = args.output

# Universe with pdb file
z = mda.Universe(reference_pdb_file)

# Coordinates are read from AVERAGE-weighted-pdb.f90 output
coordinates = read_coordinate_data(coordinate_file)

# PDB atoms
selected_atoms = z.atoms

# Replace the positions of the atoms in the pdb file with the atoms from the f90 average structure
selected_atoms.positions = coordinates

# Write out universe with correct atoms
selected_atoms.write(output_pdb_file)
