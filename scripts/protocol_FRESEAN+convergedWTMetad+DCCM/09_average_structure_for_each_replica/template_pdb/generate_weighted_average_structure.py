# Header
#  
# In this code we will generate the weighted average structure from
# a provided WT- metadynamics trajectory (07-metadyn/) and weights (08-reweight/). 
#

import MDAnalysis as mda
import numpy as np
import math
import argparse

#Parse arguments for reference pdb and the trajectory
parser = argparse.ArgumentParser()
parser.add_argument("-ref", "--reference", type=str, required=True)
parser.add_argument("-trr", "--trajectory", type=str, required=True)
parser.add_argument("-cv", "--colvar", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

reference_pdb = args.reference
trajectory = args.trajectory
colvar_file = args.colvar
output_file = args.output

# Import trajectory
u = mda.Universe(reference_pdb, trajectory)

# All atoms are used to compute average structure
ag = u.select_atoms("all")
num_frames = u.trajectory.n_frames 
nAtoms = len(ag) # Number of atoms

# Empty array used to store final structure
average_position = np.zeros((nAtoms, 3))

# Parameters
bKT=2.5 
wf=10e25 #Rescaling 
skip=10 # Compute every 10 frames


# Weights are extracted from the final column in the reweight file
wk = np.exp(np.loadtxt(colvar_file, comments="#")[::skip,-1]/bKT)/wf
wt = np.sum(wk)

# For each frame in the trajectory
for frame in u.trajectory[::skip]:
    current_frame_idx = int(frame.frame/skip) # Get the current timestep index

    # Add to the average position
    average_position = average_position + wk[current_frame_idx]*ag.positions

# Write the final weighted average structure in pdb format
with open(f"{output_file}","w") as file:
    for i in range(nAtoms):
        file.write(f"{format(average_position[i,0], '.4f')}\t{format(average_position[i,1], '.4f')}\t{format(average_position[i,2], '.4f')}\t{format(wt, '.4f')}\n")
