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
EIGVEC_D1 = np.loadtxt(f'{projDispFile}', dtype=float, comments=['#'], usecols=[1,2])

# Convert numpy array to pandas dataframe
# Convert from Angstroms to nanometers
EIGVEC_D1_df = pd.DataFrame(EIGVEC_D1, columns=["Mode 7", "Mode 8"])

range_mode7 = np.std(EIGVEC_D1_df["Mode 7"])
range_mode8 = np.std(EIGVEC_D1_df["Mode 8"])
print(f"Standard Deviation of Mode 7 Displacement Projection: {range_mode7} nm")
print(f"Standard Deviation of Mode 8 Displacement Projection: {range_mode8} nm")
