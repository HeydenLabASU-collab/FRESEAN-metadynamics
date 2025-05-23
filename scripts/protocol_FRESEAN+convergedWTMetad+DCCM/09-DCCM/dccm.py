# %%
import MDAnalysis as mda
import numpy as np
import unwrap as pbc
import align as fit
from tqdm import *

# %%
nReplicas = 20
temperatures = [300] * nReplicas

# %%
# specify input file names
TOPOL = [f"../07-metadyn/metadyn_{i}/metadyn_prot.tpr" for i in range(nReplicas)]
TRAJ = [f"../07-metadyn/metadyn_{i}/metadyn_prot_pbc.trr" for i in range(nReplicas)]
WEIGHTS = [f"../08-reweight/reweight_{i}/plumed-reweight-CV.out" for i in range(nReplicas)]
# specify output file names
AVGCRD = [f"prot_avg_replica-{i}.xyz" for i in range(nReplicas)]
AVGCRD_Ca = [f"prot-Ca_avg_replica-{i}.xyz" for i in range(nReplicas)]
AVGCRD_all = "prot_avg.xyz"
AVGCRD_Ca_all = "prot-Ca_avg.xyz"
DCCM = [f"dccm_replica-{i}.dat" for i in range(nReplicas)]
DCCM_all = "dccm_avg.dat"

# %%
# create temporary universe to fix reference coordinates
u = mda.Universe(TOPOL[0], TRAJ[0])
ca = u.select_atoms("protein and name CA")
nRes = ca.n_atoms
refPos = np.copy(ca.atoms.positions)

# %%
print(f'computing weighted average structures for each replica')
for r in range(nReplicas):
    print(f' -> replica {r+1}/{nReplicas}')
    # Load the universe
    u = mda.Universe(TOPOL[r], TRAJ[r])
    # Select protein and C-alpha atoms
    prot = u.select_atoms("protein")
    ca = u.select_atoms("protein and name CA")
    
    # read weights
    kbT = 0.0083145 * temperatures[r]
    frameWeights = np.exp(np.loadtxt(WEIGHTS[r], comments="#")[:,-1]/kbT)
    sumWeights = np.sum(frameWeights)

    # prepare unwrapping and alignment (if needed)
    unwrap = pbc.unwrap(u)
    # ensure that the reference coordinates are the same for all replicas
    # before initializing alignment
    ca.atoms.positions = refPos
    align = fit.align(u,ca)

    # array for (weighted) average positions
    avgPos = np.zeros((prot.n_atoms,3))
    frame = 0
    # loop over trajectory frames
    for ts in tqdm(u.trajectory):
        unwrap.single_frame()
        align.single_frame()
        avgPos += frameWeights[frame] * prot.positions
        frame += 1
    avgPos /= sumWeights
    # write average coordinates to file
    prot.atoms.positions = avgPos
    prot.write(AVGCRD[r], format="xyz")
    ca.write(AVGCRD_Ca[r], format="xyz")

# %%
print(f'computing average structure for all replicas')
# averaging over replicas
avgPos = np.zeros((prot.n_atoms,3))
for r in range(nReplicas):
    # Load the universe
    u = mda.Universe(AVGCRD[r], AVGCRD[r])
    avgPos += u.atoms.positions
avgPos /= nReplicas
# write average coordinates to file
u.atoms.positions = avgPos
u.atoms.write(AVGCRD_all, format="xyz")

# averaging over replicas
avgPos = np.zeros((ca.n_atoms,3))
for r in range(nReplicas):
    # Load the universe
    u = mda.Universe(AVGCRD_Ca[r], AVGCRD_Ca[r])
    avgPos += u.atoms.positions
avgPos /= nReplicas
# write average coordinates to file
u.atoms.positions = avgPos
u.atoms.write(AVGCRD_Ca_all, format="xyz")

# %%
print(f'computing weighted dynamic cross correlation matrix (DCCM) for each replica')
for r in range(nReplicas):
    print(f' -> replica {r+1}/{nReplicas}')
    # Load the universe
    u = mda.Universe(TOPOL[r], TRAJ[r])
    # Select C-alpha atoms
    ca = u.select_atoms("protein and name CA")

    # read weights
    kbT = 0.0083145 * temperatures[r]
    frameWeights = np.exp(np.loadtxt(WEIGHTS[r], comments="#")[:,-1]/kbT)
    sumWeights = np.sum(frameWeights)

    # prepare unwrapping and alignment (if needed)
    unwrap = pbc.unwrap(u)
    # ensure that the reference coordinates are the same for all replicas
    # before initializing alignment
    ca.atoms.positions = refPos
    align = fit.align(u,ca)

    # array for DCCM matrix and squared displacements
    dccm = np.zeros((nRes, nRes))
    sqDisp = np.zeros(nRes)
    frame = 0
    # loop over trajectory frames
    for ts in tqdm(u.trajectory):
        unwrap.single_frame()
        align.single_frame()
        # compute DCCM
        disp = ca.atoms.positions - avgPos
        sqDisp += frameWeights[frame] * np.sum(disp**2, axis=1)
        dccm += frameWeights[frame] * np.dot(disp, disp.T)
        frame += 1
    sqDisp /= sumWeights
    roots = np.sqrt(sqDisp)
    for i in range(nRes):
        dccm[i,:] /= roots[i]
        dccm[:,i] /= roots[i]
    dccm /= sumWeights
    # write DCCM to file
    np.savetxt(DCCM[r], dccm, fmt="%10.4f", delimiter=",")

print(f'computing average DCCM for all replicas')
avgDCCM = np.zeros((nRes, nRes))
for r in range(nReplicas):
    avgDCCM += np.loadtxt(DCCM[r], delimiter=",")
avgDCCM /= nReplicas
# write average DCCM to file
np.savetxt(DCCM_all, avgDCCM, fmt="%10.4f", delimiter=",")


