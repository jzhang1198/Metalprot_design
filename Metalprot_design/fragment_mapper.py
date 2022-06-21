"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This file contains functions for mapping designed fragments back to library fragments.
"""

#imports
import os
import pickle
import numpy as np
import pandas as pd
from prody import AtomGroup, parsePDB

def _get_neighbors(structure, coordinating_resind: int, no_neighbors: int):
    """Helper function for extract_cores. Finds neighbors of an input coordinating residue.

    Args:
        structure (prody.AtomGroup): Structure of full-length protein.
        coordinating_resind (int): Resindex of coordinating residue.
        no_neighbors (int): Number of neighbors in primary sequence to coordinating residue.

    Returns:
        core_fragment (list): List containing resindices of coordinating residue and neighbors. 
    """

    chain_id = list(set(structure.select(f'resindex {coordinating_resind}').getChids()))[0]
    all_resinds = structure.select(f'chain {chain_id}').select('protein').getResindices()
    terminal = max(all_resinds)
    start = min(all_resinds)

    extend = np.array(range(-no_neighbors, no_neighbors+1))
    _core_fragment = np.full((1,len(extend)), coordinating_resind) + extend
    core_fragment = [ind for ind in list(_core_fragment[ (_core_fragment >= start) & (_core_fragment <= terminal) ]) if ind in all_resinds] #remove nonexisting neighbor residues

    return core_fragment

def construct_library(core_database: str):
    
    for pdb in [os.path.join(core_database, file) for file in os.listdir(core_database) if 'pdb' in file]:
        core = parsePDB(pdb)
        metal_resindex = core.select('hetero').select('name NI MN ZN CO CU MG FE').getResindices()[0]
        coordinating_resindices = list(set(core.select(f'protein and not carbon and not hydrogen and within 2.83 of resindex {metal_resindex}').getResindices()))

        for resindex in coordinating_resindices:
            fragment = _get_neighbors(core, resindex, 1)
            



    pass


def _load_library(compiled_features_file: str):
    with open(compiled_features_file, 'rb') as f:
        library = pickle.load(f)

    return library

def map_fragments(core: AtomGroup, compiled_features_file: str):

    library = _load_library(compiled_features_file)
    distance_matrices = np.vstack(list(library['distance_matrices']))





