# FRESEAN-metadynamics
This codebase allows the user to run MD simulations of proteins, extract vibrational motions using the FREquency-SElective ANharmonic (FRESEAN) mode analysis, and to use the lowest-frequency vibrational modes as collective variables in enhanced sampling simulations (well-tempered metadynamics) to speed up conformational sampling. 
Please read through each of the following sections to understand how to use this repository. 
<details>
  
<summary> FRESEAN Installation </summary>

## Tested Dependencies
- FFTW 3.3.10 (https://www.fftw.org/doc/Installation-and-Customization.html)
- GNU make 3.8.1 (https://www.gnu.org/software/make/)
- gcc 14.2 (https://gcc.gnu.org/gcc-14/)
- python 3.12 (https://www.python.org/downloads/)
- gromacs 2022.5 (https://manual.gromacs.org/2022.5/download.html)
- plumed 2.8 (https://www.plumed.org/download)
- conda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

## FRESEAN Installation
Please follow the following instruction to install our suite of tools.
```
git clone https://github.com/HeydenLabASU-collab/FRESEAN-metadynamics.git
cd FRESEAN-metadynamics
make
make install
make clean
conda env create -f environment.yml
source ~/.bashrc
```
## FRESEAN Usage
- To get started, run `fresean` on the command line to see available tools. <br>
- To understand required input formats, run `fresean -h` <br>
-To see detailed information about each subroutine, run `fresean subroutine -h`, where `subroutine` is replaced with the subroutine name. <br><br>

If you have already set up GROMACS 2022.5 with Plumed 2.8, please proceed to the **FRESEAN Toolbox Programs** section to get an overview of the provided tools. Proceed to the **Provided Protocols-HEWL Example** section to get an overview on provided scripts.
</details>

<details>
  
<summary> Setting up the MD Engine </summary>


# Setting up the MD Engine

Gromacs 2022.5 is used as the MD engine and Plumed 2.8 is used as a plugin to run metadynamics.

## Step 1: Compiling PLUMED 2.8
Download PLUMED 2.8 from here: https://www.plumed.org/download
```
interactive
tar xfz plumed-2.8.tgz
cd plumed-2.8
./configure --prefix=$HOME/plumed-2.8
make -j 4
make install
```

Make sure that these paths are included in your `.bashrc` file otherwise `plumed` won't be found.
## Step 1B: BASHRC FILE FOR SOURCING
```
export PATH=$PATH:$HOME/plumed-2.8/bin
export PLUMED_VIMPATH=$PLUMED_VIMPATH:$HOME/plumed-2.8/lib/plumed/vim
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$HOME/plumed-2.8/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/plumed-2.8/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/plumed-2.8/lib/pkgconfig
export PLUMED_KERNEL="$HOME/plumed-2.8/lib/libplumedKernel.so"
```

Once plumed is installed, gromacs can be installed with the plumed patch. If patching returns `--runtime not found`, make sure that there are two dashes in front of `runtime`.
## Step 2: Patching PLUMED 2.8 in Gromacs 2022.5
Download GROMACS 2022.5 from here: https://manual.gromacs.org/documentation/2022.5/download.html
```
cd ..
tar xfz gromacs-2022.5.tar.gz
mv gromacs-2022.5 gromacs-2022.5-plumed-2.8
cd gromacs-2022.5-plumed-2.8
plumed patch -p --runtime
cd ..
```

## Step 3: Compiling Gromacs 2022.5 with PLUMED 2.8 (example for the SOL cluster at Arizona State University)
```
module load gcc-11.2.0-gcc-11.2.0
module load cuda-11.7.0-gcc-11.2.0
cd gromacs-2022.5-plumed-2.8
mkdir build
cd build
cmake .. -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2022.5-plumed-2.8 -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_plumed -DGMX_LIBS_SUFFIX=_plumed -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
make -j 4
make check
make install
```
You can now use GROMACS with plumed with the command `gmx_plumed`. I also have just normal gromacs installed, which is why I have the weird name change. If you just have gromacs with plumed, the cmake command should be modified to `cmake .. -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2022.5-plumed-2.8 -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON`.

## Configuring your BASHRC
Add the following line to your `~/.bashrc` file. Don't forget to run `source ~/.bashrc`!
```
source '$HOME/gromacs-2022.5-plumed-2.8/bin/GMXRC.bash'
```
</details>

<details>
<summary> Provided Protocols - HEWL Example </summary>

Two protocols are provided in the `scripts` directory. 

1. `protocol_FRESEAN+singleWTMetad`: Starting from a pdb file, run a single 250 ns well-tempered metadynamics with 0 THz FRESEAN modes as collective variables.
2. `protocol_FRESEAN+convergedWTMetad`: Starting from a pdb file, run 20x100 ns well-tempered metadynamics with 0 THz FRESEAN modes as collective variables. __This protocol will generate converged free energy surfaces.__

# Protocol for Single Well-Tempered Metadynamics Simulation
There are example scripts provided at `scripts/protocol_FRESEAN+singleWTMetad`. This workflow starts with a pdb file and runs 250 ns well-tempered metadynamics with 0 THz FRESEAN modes as collective variables. There is a `run.sh` script in each folder that runs the respective step. Each step is briefly described below.

## 00-prep/run.sh 
Prepare your simulation by adding a box around the protein, adding solvent, and generating ions. Keep in mind that the pdb filename, force field, water model and box size will have to be set manually. The default is force field used in our scripts is AMBER99sb-ILDN and tip3p.

## 01-em+equi/run.sh 

Energy minimization (`em.mdp`) and 100ps NPT equilibration (`equi.mdp`). 

## 02-MD/run.sh

20 ns NPT sampling simulation (`sample-NPT.mdp`) with 20 fs output frequency.

## 03-CG/run.sh

Coarse-grain simulation using `fresean coarse`.

## 04-FRESEAN/run.sh

Generate velocity cross-correlation matrix (`fresean matrix`), diagonalize the matrix (`fresean eigen`) and extract 0 THz Modes 7 and 8 (`fresean extract`) to `.xyz` format.

## 05-ModeProj/run.sh

Displacement projection of 20 ns trajectory onto FRESEAN modes. Conversion of coarse-grain FRESEAN modes to an all-atom representation (`plumed-mode-input.pdb`).

## 06-metadyn/run.sh

Run 250 ns NPT WT-metadyn simulation with plumed input file `plumed-mode-metadyn.dat`. Hills file will be `plumed-mode-metadyn.hills`. Calculate 2D free energy surface in FRESEAN space (`plumed-mode-metadyn.fes`).

## 07-reweight/run.sh

Reweight 2D free energy surface in FRESEAN space to new collective variable space. This will require a plumed input file (`07-reweight/plumed-reweight-CV.dat`) where you can define the space you are reweighting into. More information on construction of the plumed input file format can be found in the PLUMED documentation (https://www.plumed.org/doc-v2.8/user-doc/html/).

# Protocol for Converged Well-Tempered Metadynamics Simulations
There are example scripts provided at `scripts/protocol_FRESEAN+convergedWTMetad`. This workflow starts with a pdb file and runs 20x100 ns well-tempered metadynamics simulations with 0 THz FRESEAN modes as collective variables. There is a `run.sh` script in each folder that runs the respective step. The first 5 steps (`00-prep` through `05-modeProj`) are the same as the first protocol. Unique steps (`06-resample` through `08-reweight`) are described below.

## 06-resample/run.sh

100 ns NPT sampling simulation (`sample-states.mdp`). Protein configurations are extracted every 5 ns.

## 07-metadyn/run.sh

Run 20x100 ns NPT WT-metadyn simulation with plumed input file `plumed-mode-metadyn.dat`. Hills file will be `plumed-mode-metadyn.hills`. Calculate 2D free energy surface in FRESEAN space (`plumed-mode-metadyn.fes`).

## 08-reweight/run.sh

Reweight 2D free energy surface in FRESEAN space to new collective variable space. This will require a plumed input file (`08-reweight/plumed-reweight-CV.dat`) where you can define the space you are reweighting into. More information on construction of the plumed input file format can be found in the PLUMED documentation (https://www.plumed.org/doc-v2.8/user-doc/html/).

# Protocol for Dynamic Cross-Correlation Analysis

There are example scripts provided at `scripts/protocol_FRESEAN+convergedWTMetad+DCCM`. This workflow contains the same first 8 steps as `scripts/protocol_FRESEAN+convergedWTMetad`.  Starting from converged free energy surfaces, the unique two steps (`09-DCCM` to `10-reweight`) are described below.

## 09-DCCM/run.sh
Compute Dynamic Cross-Correlation Matrix (DCCM) as a weighted average over all replicas `07-metadyn/run.sh`. Then analyze the DCCM to extract distances between collectively moving domains and generate a PLUMED input file to compute free energy surfaces with proper reweighting. This part of the analysis is implemented in python and uses MDAnalysis. Run `conda env create --file=environment.yml` in the top directory of the FRESEAN-metadynamics repository to create the environment. Then modify 09-DCCM/run.sh to ensure the environment is activated.

## 10-reweight/run.sh
Reweight 2D free energy surface in FRESEAN space into space defind by distances between collectively moving domains identified in the DCCM. This will require a plumed input file generated by 09-DCCM/run.sh. More information on construction of the plumed input file format can be found in the PLUMED documentation (https://www.plumed.org/doc-v2.8/user-doc/html/).

</details>

<details>
  
<summary> FRESEAN Toolbox Programs </summary>

# FRESEAN Toolbox Programs

Information about the available tools are also accessible by running `fresean`.

## `fresean freqs`
`fresean freqs` will output the current frequency resolution. For example, the following command computes the frequency resolution of analysis performed with a correlation function length of 2 ps (# of Correlation Points = 100, Timestep = 0.02 ps).
```
fresean freqs -n 100 -t 0.02 -o freqs.txt
```

## `fresean mtop`
`fresean mtop` converts a GROMACS topology file (`.top`) into a custom topology format (`.mtop`) recognized by other `fresean` subroutines. 
```
fresean mtop -p complex.top
```

## `fresean coarse`
`fresean coarse` converts an all-atom trajectory into a coarse-grained trajectory consisting of center-of-mass sidechain and backbone beads. Example of the `coarse.inp` input file can be fould in the __FRESEAN Toolbox Input Files__ section.

> **Warning**
> Currently, `fresean coarse` is only suported for canonical amino acids.

> **Warning**
> Coarse-grained trajectory is output in `.gro` format and must be converted to `.trr` manually. This can be done with `gmx trjconv` (see example script `scripts/protocol_FRESEAN+singleWTMetad/03-CG/run.sh` for more info on how this can be done)
>

```
fresean coarse -f coarse.inp
```

## `fresean matrix`
`fresean matrix` generates a frequency dependent cross-correlation matrix in binary (`.mmat`) format. Example of the `matrix.inp` input file can be fould in the __FRESEAN Toolbox Input Files__ section.

```
fresean matrix -f matrix.inp
```

## `fresean eigen`
`fresean eigen` diagonalizes the frequency dependent cross-correlation matrix generated by `fresean matrix`. The `-n` parameter is set to the number of points in the correlation function. Produces a human-readable file of eigenvalues (`eval*.dat`)  and binary file of eigenvectors (`evec*.mmat`). Each row in the eigenvalue file corresponds to a different frequency, starting with zero frequency in the first row. Subsequent rows represent frequencies increasing by the frequency resolution.

```
fresean eigen -m matrix_fresean.mmat -n 100
```

## `fresean extract`
`fresean extract` converts the binary eigenvector file produced by `fresean eigen` into human-readable `xyz` format. Example of the `extract.inp` input file can be fould in the __FRESEAN Toolbox Input Files__ section.

```
fresean extract -f extract.inp
```

</details>
<details>
<summary> Using FRESEAN Modes as Collective Variables </summary>

## Using FRESEAN Modes as Collective Variables
To use the FRESEAN modes as collective variables in a well-tempered metadynamics simulation, the vibrational modes generated by FRESEAN must be:

1. Converted to a PLUMED-compatable file format (`.pdb`)
2. Converted to an all-atom representation

The script `scripts/protocol_FRESEAN+singleWTMetad/05-ModeProj/prep_plumed.py` performs both conversions. To convert a coarse-grained mode to an all-atom representation, the per-bead components of a coarse-grained FRESEAN mode are distributed over the atoms contributing to each bead such that the sum of atomic components equals the corresponding per-bead components.

> **Note**
> Plumed requires __one__ input file (in `.pdb` format) that contains a reference structure and all collective variables. More information can be found in the PLUMED documentation (https://www.plumed.org/doc-v2.8/user-doc/html/).

## Well-Tempered Metadynmaics Parameters

- Initial Gaussian Height: 0.1 kJ/mol
- Gaussian Width: 0.001
- Deposit Rate: 1 Gaussian per picosecond (ps)
- Bias Factor: 10
- Temperature: 300 K

</details>


<details>
<summary> FRESEAN Toolbox Input Files </summary>

# FRESEAN Toolbox Input Files
For most of the toolbox, input files (`.inp`) are utilized to gather dependencies. __Lines starting with a hash are ignored__ but can serve as a header for the variable on the subsequent line. Here is an example of the format.

```
#fnTop
complex.mtop
```

Input parameters depend on the subroutine.
> **Note**
> Parameters with a green ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) box are system-specific and should be modified depending on your protein and filenames. Parameters with a red ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) box are can be treated as static variables (recommended values are provided). 

## `fresean coarse` input format
This input file is utilized by `fresean coarse`.<br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnTop`: `.mtop` file generated by `fresean mtop` <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnCrd`: All-atom trajectory (only `.trr` format is tested extensively) <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnVel`: If `fnCrd` is `.trr` format, this should not be set. If trajectory is in `.xyz`/`.crd`/`.dcd` format, `fnCrd` is the file containing trajectory positions and `fnVel` is the file containing trajectory velocities. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `fnJob`: `.job` file used to define atom groups. All atoms are selected as default. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `grp`: Group number to read from `.job` file. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `nRead`: Number of frames to read from trajectory. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `analysisInterval`: Trajectory read-in frequency (1 = read in every frame). <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOutTraj`: Output coarse-grained trajectory in `.gro` format. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOutTopol`: Output coarse-grained topology in `.mtop` format. <br>

## `fresean matrix` input format
This file is utilized by `fresean matrix`. There are two files provided in `inp-files`. `02-gen-modes.inp` is used for all-atom analysis and `02-cgen-modes` is used if `fresean coarse` was run first (for proteins; recommended if goal is to extract collective variables for enhanced sampling). <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnTop`: `.mtop` file generated by `fresean mtop` <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnCrd`: All-atom trajectory (`.trr` recommended) <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnVel`: If `fnCrd` is `.trr` format, this should not be set. If trajectory is in `.xyz`/`.crd`/`.dcd` format, `fnCrd` is the file containing trajectory positions and `fnVel` is the file containing trajectory velocities. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `fnJob`: `.job` File used to define atom groups that can be selected later as `alignGrp` and `analyzeGrp`. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `nRead`: Number of frames to read from trajectory. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `analysisInterval`: Trajectory read-in frequency (1 = read in every frame). <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnRef`: Reference file used for translational and rotational fitting <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `alignGrp`: Index of atom group for translational and rotational alignment. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `analyzeGrp`: Index of atom group for analysis. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `wrap`: Treatment of periodic boundary conditions for distance calculations. Matters only for dynamic distance-based selections (not expected here). Set to 0 (do nothing). <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `nCorr`: Number of delay times in the correlation function. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `winSigma`: width in wavenumbers (cm^-1) of Gaussian filter applied during FFT. Recommended value is 10 cm^-1.<br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `binaryMatrix`: Format of matrix output. 0 for ASCII, 1 for binary `.mmat`. Recommend option `1` due to storage cost of ASCII format. <br>
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `doGenModes`: Compute generalized normal modes through minimization of off-diagonal matrix elements simultaenously at all frequencies as described in G. Mathias & M. Baer, J. Chem. Theory Comput. 2011, 7, 2028-2039. Set to `0` if using [FRESEAN Toolbox Workflow](#FRESEAN-Toolbox-Workflow).
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `convergence`: Convergence criteria for generalized normal mode calculation. 
- ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) `maxIter`: Maximum number of swaps for generalized normal mode calculation.
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOut`: Name of output `.mmat` file.


## `fresean extract` input format
This input file is utilized by `fresean extract`.
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnEigVec`: Binary eigenvector (`.mmat`) file. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `extractMode`: If set to option `0`, `freqSel` will extract the eigenvectors at the nearest frequency (in wavenumbers). If set to option `1`, `freqSel` will extract the eigenvectors at the provided matrix index. Mode 1 is recommended. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `freqSel`: If `extractMode` set to option `0`, this is the desired extraction __frequency (in wavenumbers)__. If `extractMode` set to 1, this is the desired extraction __matrix index__. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `trajSavingStep`: Time between time frames in trajectory (in picoseconds). <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `modeStart`: First mode to extract. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `modeEnd`: Last mode to extract. <br>
- ![#c5f015](https://placehold.co/15x15/c5f015/c5f015.png) `fnOut`: String to append to output file. <br>
</details>

<details>

<summary> References </summary>

# References
- M. A. Sauer, S. Mondal, B. Neff, S. Maiti, M. Heyden, Fast Sampling of Protein Conformational Dynamics, arXiv:2411.08154 
- S. Mondal, M. A. Sauer, M. Heyden, Exploring Conformational Landscapes Along Anharmonic Low-Frequency Vibrations, J. Phys. Chem. B 128, 7112-7120 (2024).
- M. A. Sauer, M.Heyden, Frequency-Selective Anharmonic Mode Analysis of Thermally Excited Vibrations in Proteins, J. Chem. Theory Comput. 19, 5481-5490 (2023).
