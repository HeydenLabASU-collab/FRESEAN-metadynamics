import ctypes as ct
import numpy as np
import MDAnalysis as mda
import time

class t_branch(ct.Structure):
    """recursive datatype to describe molecules as trees"""
    pass
t_branch._fields_ = (("node", ct.c_int32),
                   ("links",ct.POINTER(t_branch)),
                   ("nLinks",ct.c_int32))

class t_trees(ct.Structure):
    """array of molecule trees"""
    _fields_ = (("trees",ct.POINTER(t_branch)),
                 ("nTrees",ct.c_int32))

#load the shared library with C routines
clib = ct.cdll.LoadLibrary("ctypes/libunwrap.so")

#define argument types of function 'buildTrees' in imported library 'clib'
clib.buildTrees.argtypes = [
    ct.POINTER(t_trees),
    ct.c_int32,
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ct.c_int32,
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
]
#define return type of function 'buildTrees' in imported library 'clib'
clib.buildTrees.restype = ct.c_int32

#define argument types of function 'unwrap' in imported library 'clib'
clib.unwrap.argtypes = [
    t_trees,
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
]
#define return type of function 'unwrap' in imported library 'clib'
clib.unwrap.restype = ct.c_int32

class unwrap:
    """make molecules whole in PBC trajectories"""
    def __init__(self,u):
        self.u = u
        self.trees = self.buildTrees()
        
    def buildTrees(self):
        """build recursive bond trees that define molecules"""
        # print('building intra-molecular bond trees ...')
        start=time.process_time()
        nAtoms=len(self.u.atoms)
        atomTags=np.zeros(nAtoms,dtype=np.int32)

        bondList=self.u.bonds.indices
        bondList.sort(axis=1)
        nBonds=len(self.u.bonds)
        p=np.zeros(nBonds,dtype=np.int64)
        for i in range(nBonds):
            #we calculate for each bond a priority for sorting
            p[i]=bondList[i][0]*nAtoms+bondList[i][1]
        order=np.argsort(p)
        #now we have a sorted list of bonds
        bondList=bondList[order]

        bondTags=np.zeros(len(self.u.bonds),dtype=np.int32)

        masses=self.u.atoms.masses.astype(np.float32)

        trees = t_trees()
        error = clib.buildTrees(
            ct.pointer(trees),
            ct.c_int(nAtoms),
            atomTags,
            ct.c_int(nBonds),
            bondTags,
            bondList,
            masses
        )
        if error != 0:
            print(f'ERROR reported by \'buildTrees\' function\nsee \'error.log\'\n')
        stop=time.process_time()
        # print(f'unwrap -> detected {trees.nTrees} molecules')
        # print(f'unwrap -> tree building time: {stop-start:.2f}s')
        # print(f'unwrap -> ready to unwrap')
        return trees

    def single_frame(self):
        """
        unwrap coordinates: 
        ensure that all components of covalent bonds are shorter than 1/2 the box
        """
        error = clib.unwrap(
            self.trees,
            self.u.trajectory.ts._pos,
            self.u.dimensions
        )
        if error != 0:
            print(f'ERROR reported by \'unwrap\' function\nsee \'error.log\'\n')
