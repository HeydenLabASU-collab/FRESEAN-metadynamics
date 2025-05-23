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

# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
dccm = np.loadtxt(DCCM_all, delimiter=",")

# %%
plt.figure(figsize=(8, 6))
im = plt.imshow(dccm, cmap='bwr', vmin=-1, vmax=1)
plt.gca().invert_yaxis()
plt.colorbar(im, label='Correlation')
plt.title('Dynamic Cross Correlation Matrix (DCCM)')
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.tight_layout()
plt.savefig("dccm.png", dpi=300)

# %%
pair1 = np.unravel_index(np.argmin(dccm), dccm.shape)
col0 = np.abs(dccm[:, pair1[0]])
col1 = np.abs(dccm[:, pair1[1]])
least_correlated = np.argmin(col0 + col1)
anti_correlated = np.argmin(dccm[:, least_correlated])
pair2=[least_correlated, anti_correlated]

# %%
domain1 = np.where(dccm[pair1[0]] > 0.5)[0]
domain2 = np.where(dccm[pair1[1]] > 0.5)[0]
domain3 = np.where(dccm[pair2[0]] > 0.5)[0]
domain4 = np.where(dccm[pair1[1]] > 0.5)[0]

# %%
plt.figure(figsize=(8,5))
x = np.arange(dccm.shape[1])
y1 = dccm[pair1[0]]
y2 = dccm[pair1[1]]

plt.plot(x, y1, marker='o', markersize=3, color='blue', label=f'DCCM row {pair1[0]+1}')
plt.plot(x, y2, marker='o', markersize=3, color='red', label=f'DCCM row {pair1[1]+1}')

# Fill area above y=0.5 for y1
plt.fill_between(x, 0.5, y1, where=(y1 > 0.5), color='blue', alpha=0.3)
# Fill area above y=0.5 for y2
plt.fill_between(x, 0.5, y2, where=(y2 > 0.5), color='red', alpha=0.3)

plt.xlabel('Residue Index')
plt.ylabel('Correlation')
plt.title('row-wise correlation')
plt.grid(True)
plt.legend(loc='center', bbox_to_anchor=(0.7, 0.1))
plt.xlim(x.min(), x.max())
plt.ylim(-0.75, 1.0)
plt.tight_layout()
plt.savefig("pair1.png", dpi=300)

# %%
plt.figure(figsize=(8,5))
x = np.arange(dccm.shape[1])
y1 = dccm[pair2[0]]
y2 = dccm[pair2[1]]

plt.plot(x, y1, marker='o', markersize=3, color='blue', label=f'DCCM row {pair2[0]+1}')
plt.plot(x, y2, marker='o', markersize=3, color='red', label=f'DCCM row {pair2[1]+1}')

# Fill area above y=0.5 for y1
plt.fill_between(x, 0.5, y1, where=(y1 > 0.5), color='blue', alpha=0.3)
# Fill area above y=0.5 for y2
plt.fill_between(x, 0.5, y2, where=(y2 > 0.5), color='red', alpha=0.3)

plt.xlabel('Residue Index')
plt.ylabel('Correlation')
plt.title('row-wise correlation')
plt.grid(True)
plt.legend(loc='center', bbox_to_anchor=(0.7, 0.1))
plt.xlim(x.min(), x.max())
plt.ylim(-0.75, 1.0)
plt.tight_layout()
plt.savefig("pair2.png", dpi=300)

# %%
def generate_plumed_distance_input(domain1, domain2, domain3, domain4):
    # Convert to 1-based indexing for PLUMED
    group1_indices = ",".join(str(i+1) for i in domain1)
    group2_indices = ",".join(str(i+1) for i in domain2)
    group3_indices = ",".join(str(i+1) for i in domain3)
    group4_indices = ",".join(str(i+1) for i in domain4)
    plumed_input = (
        "#Reference file containing eigenvectors\n"
        "PCAVARS REFERENCE=plumed-mode-input.pdb TYPE=OPTIMAL LABEL=pca\n"
        "\n"
        "# Go through metadynamics trajectory and get the weight every frame\n"
        "METAD ...\n"
        "LABEL=metad\n"
        "ARG=pca.eig-1,pca.eig-2\n"
        "PACE=10000\n"
        "HEIGHT=0.0\n"
        "BIASFACTOR=10\n"
        "SIGMA=0.001,0.001\n"
        "FILE=plumed-mode-metadyn.hills\n"
        "TEMP=300.0\n"
        "RESTART=YES\n"
        "... METAD\n"
        "\n"
        f"g1 : COM ATOMS={group1_indices}\n"
        f"g2 : COM ATOMS={group2_indices}\n"
        "d1: DISTANCE ATOMS=g1,g2\n"
        "\n"
        f"g3 : COM ATOMS={group3_indices}\n"
        f"g4 : COM ATOMS={group4_indices}\n"
        "d2: DISTANCE ATOMS=g3,g4\n"
        "\n"
        "# Print each of the desired quantities to a file. Can manually reweight.\n"
        "# Weight of frame i is given by w_i proportional to exp(V/kT). V is given in the reweight file as metad.bias.\n"
        "PRINT ARG=pca.eig-1,pca.eig-2,d1,d2,pca.residual,metad.bias FILE=plumed-reweight-CV.out STRIDE=1\n"
    )
    with open("plumed-reweight-CV.dat", "w") as f:
        f.write(plumed_input)

# %%
generate_plumed_distance_input(domain1, domain2, domain3, domain4)


