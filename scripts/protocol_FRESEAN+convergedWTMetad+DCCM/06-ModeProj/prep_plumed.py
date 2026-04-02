#!/usr/bin/python3

import numpy as np
import sys
import os 

def extract_atom_masses(mtop_path):
    """
    Extract (atom_name, mass) pairs from a GROMACS .mtop file.
    Looks between [atoms] and [sites] sections.
    """
    atoms_section = False
    result = []

    with open(mtop_path) as f:
        for line in f:
            stripped = line.strip()

            # Start of atoms section
            if stripped.startswith("[ atoms ]") or stripped.startswith("[atoms]"):
                atoms_section = True
                continue
            
            # End of atoms section
            if stripped.startswith("[ sites ]") or stripped.startswith("[sites]"):
                atoms_section = False

            # If inside atoms section, parse atom + mass
            if atoms_section:
                parts = stripped.split()
                if len(parts) >= 2:
                    name = parts[0]
                    try:
                        mass = float(parts[1])
                    except ValueError:
                        continue  # skip non-numeric mass rows
                    result.append(mass)

    return result

# Working directory where the reference pdb file can be found. Follow above steps.
# This reference structure is already in pdb format.


# Mode to convert to .pdb format
modefile = str(sys.argv[1])
topology = str(sys.argv[2])
modenum = int(sys.argv[3])
structfile = str(sys.argv[4])
if os.path.isfile(modefile):
    print("File exists!")
else:
    print("File does NOT exist.")


struct = np.loadtxt(f'{structfile}', dtype=str, comments=['END', 'REMARK', 'TITLE', 'MODEL', 'TER', 'CRYST1'], usecols=[0,1,2,3,4,5,6,7])
def read_XYZ_modefile(modefile):
    nModes = 0
    nLines = 0
    with open(modefile) as f:
        for line in f:
            if line.startswith("X"):
                nLines += 1
            if line.startswith("Frequency"):
                nModes += 1

    nAtoms = int(nLines/nModes)
    print("Number of Modes in XYZ File:", nModes)
    print("Number of Atoms per mode in XYZ File:", nAtoms)
    return nModes, nAtoms


nModes, nAtoms = read_XYZ_modefile(modefile)
eigvec = np.zeros((nModes,nAtoms,3))


with open(f'{modefile}','r') as file:
    currAtom = 0
    currMode = -1
    for line in file:
        if line.startswith("Frequency"):
            currMode += 1
            currAtom = 0
        if line.startswith("X"):
            eigvec[currMode,currAtom,:]=line.split()[1:4]
            currAtom += 1
        


masses = extract_atom_masses(topology)
masses = np.reshape(masses, (np.shape(masses)[0],1))
mass_weighted = np.empty_like(eigvec)
for i in range(eigvec.shape[0]):
    mass_weighted[i] = 100 * eigvec[i] * np.sqrt(masses)
print(np.shape(mass_weighted))
with open(f"evec_{modenum}_aa_scaled.pdb","w") as writer:
    writer.write("REMARK TYPE=DIRECTION\n")
    counter = 0
    for line in struct:
        writer.write(f'{line[0]}  {line[1]:>5} {line[2]:>4} {line[3]:1} {line[4]:>5}    {float(mass_weighted[modenum-1,counter,0]):8.3f}{float(mass_weighted[modenum-1,counter,1]):>8.3f}{float(mass_weighted[modenum-1,counter,2]):>8.3f}  1.00  0.00\n')
        counter += 1
    writer.write("END\n")
writer.close()
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
