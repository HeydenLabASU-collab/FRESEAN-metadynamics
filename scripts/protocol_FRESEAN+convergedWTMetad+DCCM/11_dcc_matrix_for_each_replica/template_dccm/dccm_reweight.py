import MDAnalysis as mda
import numpy as np
import math
import argparse

#Parse arguments for reference pdb and the trajectory
parser = argparse.ArgumentParser()
parser.add_argument("-ref", "--reference", type=str, required=True)
parser.add_argument("-trr", "--trajectory", type=str, required=True)
parser.add_argument("-cv", "--colvar", type=str, required=True)
parser.add_argument("-av", "--average_structure", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

reference_pdb = args.reference
trajectory = args.trajectory
colvar_file = args.colvar
avg_struct = args.average_structure
output_file = args.output


trajectory_universe = mda.Universe(reference_pdb, trajectory)

ag = trajectory_universe.select_atoms("name CA")
ca = trajectory_universe.select_atoms('name CA')
num_frames = trajectory_universe.trajectory.n_frames
n_ca= len(ca)

# Get the average structure
average_universe = mda.Universe(avg_struct)
alpha_carbons = average_universe.select_atoms("name CA")
av = alpha_carbons.positions



bKT=2.5
wf=10e15
skip=10
nSteps=100001

wk = np.exp(np.loadtxt(colvar_file, comments="#")[::skip, -1]/bKT)/wf
wt = np.sum(wk)

cor_xy=np.zeros((n_ca, n_ca))
cor_x=np.zeros((n_ca, n_ca))
cor_y=np.zeros((n_ca, n_ca))

coords=np.zeros((num_frames, n_ca, 3))

for frame in trajectory_universe.trajectory[::skip]:
    current_frame_idx = int(frame.frame/skip)
    coords[current_frame_idx] = ag.positions
    for residue_x in range(n_ca - 1):
        for residue_y in range(residue_x+1, n_ca):

            di = coords[current_frame_idx, residue_x] - av[residue_x]
            dj = coords[current_frame_idx, residue_y] - av[residue_y]

            cor_xy[residue_x, residue_y] = cor_xy[residue_x, residue_y] + wk[current_frame_idx] * (di[0]*dj[0] + di[1]*dj[1]+ di[2]*dj[2])
            cor_xy[residue_y, residue_x] = cor_xy[residue_x, residue_y]

            cor_x[residue_x, residue_y] = cor_x[residue_x, residue_y] + wk[current_frame_idx] * (di[0]**2 + di[1]**2 + di[2]**2)
            cor_x[residue_y, residue_x] = cor_x[residue_x, residue_y]

            cor_y[residue_x, residue_y] = cor_y[residue_x, residue_y] + wk[current_frame_idx] * (dj[0]**2 + dj[1]**2 + dj[2]**2)
            cor_y[residue_y, residue_x] = cor_y[residue_x, residue_y]
           

with open(output_file,"w") as file:
    for residue_x in range(n_ca):
        for residue_y in range(n_ca):
            if residue_x == residue_y: 
                cor_xy[residue_x, residue_y]=1.0
            else:
                cor_xy[residue_x, residue_y]=cor_xy[residue_x, residue_y]/((np.sqrt(cor_x[residue_x, residue_y]))*(np.sqrt(cor_y[residue_x, residue_y])))
            file.write(f"{residue_x+1}\t{residue_y+1}\t{cor_xy[residue_x, residue_y]}\n")
        file.write("\n")
