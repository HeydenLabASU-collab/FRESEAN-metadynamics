import ctypes as ct
import MDAnalysis as mda
from MDAnalysis.coordinates import core
import numpy as np

# %%
class t_residue(ct.Structure):
    '''datatype for group of atoms'''
    _fields_ = (("offset", ct.c_int32),
                ("nAtoms", ct.c_int32),
                ("masses", ct.POINTER(ct.c_double)),
                ("totMass", ct.c_double),
                ("COM", ct.c_float*3)
                )

class t_align(ct.Structure):
    '''datatype with settings for trans+rot alignment'''
    _fields_ = (("align", ct.c_int32),
                ("centerCOM", ct.c_int32),
                ("subCOMvel", ct.c_int32),
                ("placeCOMInBox", ct.c_int32),
                ("rotate", ct.c_int32),
                ("rotVel", ct.c_int32),
                ("nAtomsRef", ct.c_int32),
                ("massesRef", ct.POINTER(ct.c_double)),
                ("totMassRef", ct.c_double),
                ("crdRefFix", ct.POINTER(ct.c_float)),
                ("COMrefFix", ct.c_float*3),
                ("crdRef", ct.POINTER(ct.c_float)),
                ("nAtomsSys", ct.c_int32),
                ("crdSys", ct.POINTER(ct.c_float)),
                ("velSys", ct.POINTER(ct.c_float)),
                ("nResSys", ct.c_int32),
                ("resSys", ct.POINTER(t_residue)),
                ("box", ct.c_float*9)
                )

# %%
#load the shared library with C routines
alignLib = ct.cdll.LoadLibrary("ctypes/libalign.so")

# %%
alignLib.setup.argtypes = [
    # alignment settings
    ct.POINTER(t_align),
    # number of atoms in reference group
    ct.c_int32,
    # masses of atoms in reference group
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    # referece coordinates
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    # number of atoms in system
    ct.c_int32,
    # masses of atoms in system
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    # number of residues in system
    ct.c_int32,
    # array with number of atoms per residue in system
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    # box vectors
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
]
alignLib.setup.restype = ct.c_int32

alignLib.alignCrd.argtypes = [
    # alignment settings
    ct.POINTER(t_align),
    # coordinates of reference group
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    # coordinates of system
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    # box vectors
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
]
alignLib.alignCrd.restype = ct.c_int32

alignLib.alignCrdVel.argtypes = [
    # alignment settings
    ct.POINTER(t_align),
    # coordinates of reference group
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    # velocities of reference group
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    # coordinates of system
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    # velocities of system
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    # box vectors
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
]
alignLib.alignCrdVel.restype = ct.c_int32

# %%
class align:
    def __init__(self,
                 u,
                 refSel,
                 align=1,
                 centerCOM=1,
                 subCOMvel=0,
                 placeCOMInBox=1,
                 rotate=1,
                 rotVel=0
                 ):
        self.u = u
        self.refSel = refSel
        self.align = t_align()
        self.align.align = align
        self.align.centerCOM = centerCOM
        self.align.subCOMvel = subCOMvel
        self.align.placeCOMInBox = placeCOMInBox
        self.align.rotate = rotate
        self.align.rotVel = rotVel
        self.box = core.triclinic_vectors(self.u.dimensions)
        self.prep()

    def prep(self):
        nAtomsPerRes = np.zeros(len(self.u.residues), dtype=np.int32)
        for i, res in enumerate(self.u.residues):
            nAtomsPerRes[i] = len(res.atoms)
        error = alignLib.setup(ct.byref(self.align),
                       self.refSel.n_atoms,
                       self.refSel.atoms.masses,
                       self.refSel.atoms.positions,
                       len(self.u.atoms),
                       self.u.atoms.masses,
                       len(self.u.residues),
                       nAtomsPerRes,
                       core.triclinic_vectors(self.u.dimensions))

    def single_frame(self):
        if self.align.subCOMvel ==1 or self.align.rotVel == 1:
            alignLib.alignCrdVel(ct.byref(self.align),
                        self.refSel.atoms.positions,
                        self.refSel.atoms.velocities,
                        self.u.trajectory.ts._pos,
                        self.u.trajectory.ts._velocities,
                        core.triclinic_vectors(self.u.trajectory.ts._unitcell))
        else:
            alignLib.alignCrd(ct.byref(self.align),
                        self.refSel.atoms.positions,
                        self.u.trajectory.ts._pos,
                        core.triclinic_vectors(self.u.trajectory.ts._unitcell))
        self.box = np.array(self.align.box[0:9]).reshape(3,3)