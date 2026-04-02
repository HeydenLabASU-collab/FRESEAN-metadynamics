#!/usr/bin/python3

import numpy as np
import sys
# Working directory where the reference pdb file can be found. Follow above steps.
# This reference structure is already in pdb format.

struct = np.loadtxt(f'ref.pdb', dtype=str, comments=['END', 'REMARK', 'TITLE', 'MODEL', 'TER', 'CRYST1'], usecols=[0,1,2,3,4,5,6,7])
mode_writer = np.zeros((np.shape(struct)[0], 3))

# Working directory where the modes are calculated from

# Mode to convert to .pdb format
mode = int(sys.argv[1])

# Number of coarsened beads and array for eigenvector
nBeads = 0

# Index for selected eigenvector line (DO NOT CHANGE - SHOULD BE 0)
curr = 0

# Open the xyz eigenvector file and iterate through, storing the eigenvector for each
# atom in the array
with open(f'input-cg-modes.xyz','r') as file:

    # Index for file line
    counter = 0
    for line in file:
        if counter == 0:
            nBeads = int(line)
            eigvec = np.zeros((nBeads,3))
        counter = counter + 1

        # Skip the two line header and record the 1st through 3rd columns until you reach
        # The next mode
        if(counter > (nBeads+2)*(mode - 1)+2 and counter <= (nBeads+2)*(mode)):
            eigvec[curr,:]=line.split()[1:4]
            curr = curr + 1



print(f"Number of beads: {nBeads}\nConverting Mode {mode}\n")
# Some functions that we use to check if the current atom in the reference structure
# 1. Is part of a residue only containing backbone atoms
# 2. Is currently a backbone atom

# This will have to be expanded to consider things like capping residues,
# non-canonical AA's, etc.
# But none of the current proteins we are analyzing have those.

backboneAtoms = ["CA", "OC1", "OC2", "C", "O", "HA", "N", "H", "H1","H2","H3"]

atom_mass = {
    "H": 1.0080,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "S": 32.070,
    "P": 30.97376200
}

def atom_masses(strings):
    masses = []
    for s in strings:
        atom = s[0]          # first character of the string
        mass = atom_mass.get(atom)
        masses.append(mass)
    return masses
    
def isBackboneOnly(resName):
    if(resName == "GLY"):
        return True
    else:
        return False

def inBackbone(atomName):
    if atomName in backboneAtoms:
        return True
    else:
        return False

def calcSizesandMasses(struct, mass, nBeads):
    prevRes = int(struct[0,4]) # Resid of the first residue in pdb file
    curr = 0
    scaling = 100 # Constant to scale modes by
    beadCounter = np.zeros(nBeads)
    massCounter = np.zeros(nBeads)
    for i in range(len(struct)):

        # Check to see if the current atom is in the same residue as the last atom
        if(int(struct[i,4]) == prevRes):
            # If it is, check to see if it is an atom with only backbone atoms
            if(isBackboneOnly(struct[i,3])):
                beadCounter[curr] = beadCounter[curr] + 1
                massCounter[curr] = massCounter[curr] + mass[i]
            elif(not isBackboneOnly(struct[i,3])):
                if(inBackbone(struct[i,2])):
                    beadCounter[curr] = beadCounter[curr] + 1
                    massCounter[curr] = massCounter[curr] + mass[i]
                elif(not inBackbone(struct[i,2])):
                    beadCounter[curr+1] = beadCounter[curr+1] + 1
                    massCounter[curr+1] = massCounter[curr+1] + mass[i]

        elif(int(struct[i,4]) != prevRes):
            prevRes = int(struct[i,4])
        
            # If the last residue only contained backbone atoms, then we only need to iterate
            # once. However, the the last residue had BACK and SIDE reßßßsidues, we need to iterate
            # twice to move to the eigenvectors associated with the next residue.
            if(isBackboneOnly(struct[i-1,3])):
                curr = curr + 1
            else:
                curr = curr + 2
            
            # This is the same as the outermost if statement. This is required to handle the edge
            # case of switching residues.
            # Could probably be in a function.
            if(isBackboneOnly(struct[i,3])):
                beadCounter[curr] = beadCounter[curr] + 1
                massCounter[curr] = massCounter[curr] + mass[i]
            else:
                if(int(struct[i,4]) == prevRes and inBackbone(struct[i,2])):
                    beadCounter[curr] = beadCounter[curr] + 1
                    massCounter[curr] = massCounter[curr] + mass[i]
                elif(int(struct[i,4]) == prevRes and not inBackbone(struct[i,2])):
                    beadCounter[curr+1] = beadCounter[curr+1] + 1
                    massCounter[curr+1] = massCounter[curr+1] + mass[i]
    return beadCounter,massCounter

# Get atom masses
pdb_atomnames = struct[:,2]
pdb_atommasses = atom_masses(pdb_atomnames)
sqrtpdb_atommasses = np.sqrt(np.array(pdb_atommasses))
# Now we will parse the reference structure and match our bead eigenvectors to the atoms


prevRes = int(struct[0,4]) # Resid of the first residue in pdb file
curr = 0
scaling = 100 # Constant to scale modes by
beadCounter, massCounter = calcSizesandMasses(struct, pdb_atommasses, nBeads)
print(beadCounter, len(beadCounter))

# For all atoms in the protein
for i in range(len(struct)):

    # Check to see if the current atom is in the same residue as the last atom
    if(int(struct[i,4]) == prevRes):
        
        # If it is, check to see if it is an atom with only backbone atoms
        if(isBackboneOnly(struct[i,3])):

            # If it is backbone only, apply current coarsened eigenvector to the current atom
            mode_writer[i,0] = float(sqrtpdb_atommasses[i]*eigvec[curr,0]/np.sqrt(massCounter[curr]))
            mode_writer[i,1] = float(sqrtpdb_atommasses[i]*eigvec[curr,1]/np.sqrt(massCounter[curr]))
            mode_writer[i,2] = float(sqrtpdb_atommasses[i]*eigvec[curr,2]/np.sqrt(massCounter[curr]))

        # If it is NOT backbone only, we need to figure out if the current atom
        # is a backbone atom or a sidechain atom
        elif(not isBackboneOnly(struct[i,3])):
            
            # If it is a backbone atom, apply the current coarsened eigenvector to 
            # current atom
            if(inBackbone(struct[i,2])):
                mode_writer[i,0] = float(sqrtpdb_atommasses[i]*eigvec[curr,0]/np.sqrt(massCounter[curr]))
                mode_writer[i,1] = float(sqrtpdb_atommasses[i]*eigvec[curr,1]/np.sqrt(massCounter[curr]))
                mode_writer[i,2] = float(sqrtpdb_atommasses[i]*eigvec[curr,2]/np.sqrt(massCounter[curr]))
                
            # If it is a sidechain atom, we have to look at the next coarsened eigenvector
            # Since our convention is BACK then SIDE
            elif(not inBackbone(struct[i,2])):
                mode_writer[i,0] = float(sqrtpdb_atommasses[i]*eigvec[curr+1,0]/np.sqrt(massCounter[curr+1]))
                mode_writer[i,1] = float(sqrtpdb_atommasses[i]*eigvec[curr+1,1]/np.sqrt(massCounter[curr+1]))
                mode_writer[i,2] = float(sqrtpdb_atommasses[i]*eigvec[curr+1,2]/np.sqrt(massCounter[curr+1]))
            
    # If the current atom residue number is not equal to the residue number of the 
    # previous atom, we need to get the eigenvectors for the next residue
    elif(int(struct[i,4]) != prevRes):
        
        # The previous residue is going to be the residue of the current atom in the next
        # iteration
        prevRes = int(struct[i,4])
        
        # If the last residue only contained backbone atoms, then we only need to iterate
        # once. However, the the last residue had BACK and SIDE residues, we need to iterate
        # twice to move to the eigenvectors associated with the next residue.
        if(isBackboneOnly(struct[i-1,3])):
            curr = curr + 1
        else:
            curr = curr + 2
            
        # This is the same as the outermost if statement. This is required to handle the edge
        # case of switching residues.
        # Could probably be in a function.
        if(isBackboneOnly(struct[i,3])):
            mode_writer[i,0] = float(sqrtpdb_atommasses[i]*eigvec[curr,0]/np.sqrt(massCounter[curr]))
            mode_writer[i,1] = float(sqrtpdb_atommasses[i]*eigvec[curr,1]/np.sqrt(massCounter[curr]))
            mode_writer[i,2] = float(sqrtpdb_atommasses[i]*eigvec[curr,2]/np.sqrt(massCounter[curr]))
        else:
            if(int(struct[i,4]) == prevRes and inBackbone(struct[i,2])):
                mode_writer[i,0] = float(sqrtpdb_atommasses[i]*eigvec[curr,0]/np.sqrt(massCounter[curr]))
                mode_writer[i,1] = float(sqrtpdb_atommasses[i]*eigvec[curr,1]/np.sqrt(massCounter[curr]))
                mode_writer[i,2] = float(sqrtpdb_atommasses[i]*eigvec[curr,2]/np.sqrt(massCounter[curr]))
            elif(int(struct[i,4]) == prevRes and not inBackbone(struct[i,2])):
                beadCounter[curr+1] = beadCounter[curr+1] + 1
                mode_writer[i,0] = float(sqrtpdb_atommasses[i]*eigvec[curr+1,0]/np.sqrt(massCounter[curr+1]))
                mode_writer[i,1] = float(sqrtpdb_atommasses[i]*eigvec[curr+1,1]/np.sqrt(massCounter[curr+1]))
                mode_writer[i,2] = float(sqrtpdb_atommasses[i]*eigvec[curr+1,2]/np.sqrt(massCounter[curr+1]))


# Export the eigenvector
with open(f"evec_{mode}_aa_backmapped.xyz","w") as writer:
    writer.write(f"{nBeads}\n")
    writer.write(f"Frequency\n")
    for line in mode_writer:
        writer.write(f'X {float(line[0]):.6e} {float(line[1]):.6e} {float(line[2]):.6e}\n')
writer.close()
