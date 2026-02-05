#!/usr/bin/python3

import sys

# Variable to be set by the user
import numpy as np
import pandas as pd

num_modes = 2 # Number of modes to read in


projDispFile = str(sys.argv[1])
pts=int(sys.argv[2])

# Read in data and store in dataframe
EIGVEC_D1 = np.zeros((num_modes,pts))
EIGVEC_D1 = np.loadtxt(projDispFile, dtype=float, comments=['#'], usecols=[1,2])

# Convert numpy array to pandas dataframe
# Convert from Angstroms to nanometers
EIGVEC_D1_df = pd.DataFrame(EIGVEC_D1, columns=["Mode 7", "Mode 8"])

range_mode7 = np.std(EIGVEC_D1_df["Mode 7"])/2
range_mode8 = np.std(EIGVEC_D1_df["Mode 8"])/2
print(f'Half of Standard Deviation of Mode 7 Displacement Projection (To be used for WT-MetaD sigma): {range_mode7:.9f} nm')
print(f'Half of Standard Deviation of Mode 8 Displacement Projection (To be used for WT-MetaD sigma): {range_mode8:.9f} nm')

# writes two lines to a text file
with open("std_dev.txt", "w") as f:
    f.write(f"{range_mode7:.9f}\n")
    f.write(f"{range_mode8:.9f}\n")
