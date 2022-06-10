"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This file contains functions for loading structures for metal binding predictions.
"""

#imports 
from pypivoter.degeneracy_cliques import enumerateCliques, countCliques
from prody import parsePDB, buildDistMatrix, AtomGroup
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import numpy as np
import pandas as pd
import itertools
import operator

#imports 
from pypivoter.degeneracy_cliques import enumerateCliques
from prody import parsePDB, buildDistMatrix, AtomGroup
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import numpy as np
import pandas as pd
import itertools
import operator

def _get_helices(pdb_file: str, c_alphas: AtomGroup):
    dssp_dict, _ = dssp_dict_from_pdb_file(pdb_file)
    helical_residues = [(resnum, chid) for resnum, chid in zip(c_alphas.getResnums(), c_alphas.getChids()) if dssp_dict[(chid, (' ', int(resnum), ' '))][1] in ['H', 'G', 'I']]
    helical_residues_by_chid = [list(group) for key, group in itertools.groupby(helical_residues, operator.itemgetter(1))]

    contiguous_helices = []
    for resis in helical_residues_by_chid:
        chid = resis[0][1]
        resnums = np.array([tup[0] for tup in resis])
        discontinuities = np.concatenate((np.array([0]), np.add(np.argwhere(resnums[1:] - resnums[0:-1] != 1).flatten(), 1), np.array([resnums[-1]])))
        contiguous_helices += [[(num, chid) for num in resnums[i:j]] for i,j in zip(discontinuities[0:-1], discontinuities[1:])]
    
    dict = {}
    for ind in range(0, len(contiguous_helices)):
        for identifier in contiguous_helices[ind]:
            dict[identifier] = ind

    return dict

def _get_neighbors(structure, resind: int, no_neighbors: int):
    """Helper function for extract_cores. Finds neighbors of an input coordinating residue.

    Args:
        structure (prody.AtomGroup): Structure of full-length protein.
        coordinating_resind (int): Resindex of reference residue.
        no_neighbors (int): Number of neighbors in primary sequence to coordinating residue.

    Returns:
        core_fragment (list): List containing resindices of coordinating residue and neighbors. 
    """

    chain_id = list(set(structure.select(f'resindex {resind}').getChids()))[0]
    all_resinds = structure.select(f'chain {chain_id}').select('protein').getResindices()
    terminal = max(all_resinds)
    start = min(all_resinds)

    extend = np.array(range(-no_neighbors, no_neighbors+1))
    _fragment = np.full((1,len(extend)), resind) + extend
    fragment = [ind for ind in list(_fragment[ (_fragment >= start) & (_fragment <= terminal) ]) if ind in all_resinds] #remove nonexisting neighbor residues

    return fragment

def _compute_ca_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray):
    ca_n = np.vstack([(1 / np.linalg.norm(n - ca, axis=1))]*3).T * (n - ca)
    ca_c = np.vstack([(1 / np.linalg.norm(c - ca, axis=1))]*3).T * (c - ca)
    n1 = (ca_n + ca_c) * -1
    n2 = np.cross(ca_n, ca_c, axis=1)
    n1 = np.vstack([(1 / np.linalg.norm(n1, axis=1))]*3).T * n1
    n2 = np.vstack([(1 / np.linalg.norm(n2, axis=1))]*3).T * n2
    d = (1.54*np.sin(np.deg2rad(54.75))) * n2
    v = (1.54*np.cos(np.deg2rad(54.75))) * n1
    ca_cb = d+v

    return ca_cb

def _compute_angles(vec1: np.ndarray, vec2: np.ndarray):
    dot = np.sum(vec1 * vec2, axis=1)
    norm1, norm2 = np.linalg.norm(vec1, axis=1), np.linalg.norm(vec2, axis=1)
    angles = np.degrees(np.arccos(dot / (norm1 * norm2)))
    return angles

def _filter(adjacency_list: np.ndarray, structure: AtomGroup, distances: np.ndarray):
    all_resindices = set(np.concatenate(list(adjacency_list)))
    coordinates = dict([(resindex, structure.select('protein').select('name C CA N').select(f'resindex {resindex}').getCoords()) for resindex in all_resindices])

    n_i, n_j = np.vstack([coordinates[resindex][0].flatten() for resindex in adjacency_list[:,0]]), np.vstack([coordinates[resindex][0].flatten() for resindex in adjacency_list[:,1]])
    ca_i, ca_j = np.vstack([coordinates[resindex][1].flatten() for resindex in adjacency_list[:,0]]), np.vstack([coordinates[resindex][1].flatten() for resindex in adjacency_list[:,1]])
    c_i, c_j = np.vstack([coordinates[resindex][2].flatten() for resindex in adjacency_list[:,0]]), np.vstack([coordinates[resindex][2].flatten() for resindex in adjacency_list[:,1]])
    ca_cb_i, ca_cb_j = _compute_ca_cb(n_i, ca_i, c_i), _compute_ca_cb(n_j, ca_j, c_j)
    ca_i_ca_j = ca_j - ca_i
    angles_i, angles_j = _compute_angles(ca_cb_i, ca_i_ca_j), _compute_angles(ca_cb_j, ca_i_ca_j)

    accepted = np.argwhere(distances <= 7)
    filtered_inds = np.intersect1d(np.intersect1d(np.argwhere(angles_i < 140), np.argwhere(angles_j > 35)), np.argwhere(distances > 7))
    filtered = adjacency_list[np.union1d(accepted, filtered_inds)]
    return filtered

def identify_sites_rational(pdb_file: str, cuttoff: float, coordination_number=(2,4), no_neighbors=1):
    features, identifiers, sources = [], [], []

    max_atoms = 4 * (coordination_number[1] + (coordination_number[1] * no_neighbors * 2))
    structure = parsePDB(pdb_file)
    c_alphas = structure.select('protein').select('name CA')

    #enumerate all helical residues in the input structure
    contiguous_helices = _get_helices(pdb_file, c_alphas)
    helix_resindices = [structure.select(f'chid {tup[1]}').select(f'resnum {tup[0]}').getResindices()[0] for tup in contiguous_helices.keys()]
    helix_resindices.sort()
    helix_resindices_selstr = 'resindex ' + ' '.join([str(i) for i in helix_resindices])
    resind2id = dict([(resindex, (structure.select(f'resindex {resindex}').getResnums()[0], structure.select(f'resindex {resindex}').getChids()[0])) for resindex in set(structure.select('protein').getResindices())])

    #build alpha carbon distance matrix. iterate through columns and enumerate pairs of residues that are within a cutoff distance.
    ca_dist_mat = buildDistMatrix(c_alphas.select(helix_resindices_selstr), c_alphas.select(helix_resindices_selstr))
    adjacency_list = []
    edge_weights = np.array([])
    row_indexer = 0
    for col_ind in range(len(ca_dist_mat)):
        for row_ind in range(1+row_indexer, len(ca_dist_mat)):
            distance = ca_dist_mat[row_ind, col_ind]

            if distance <= cuttoff:
                adjacency_list.append(np.array([helix_resindices[col_ind], helix_resindices[row_ind]]))
                edge_weights = np.append(edge_weights, distance)

        row_indexer += 1
    adjacency_list = _filter(np.vstack(adjacency_list), structure, edge_weights)
    resind2map = dict([(resind, index) for index, resind in enumerate(np.sort(np.unique(adjacency_list.flatten())))])
    map2resind = dict([(index, resind) for index, resind in enumerate(np.sort(np.unique(adjacency_list.flatten())))])

    mapped_adjacency_list = np.vectorize(resind2map.get)(adjacency_list)
    cliques = enumerateCliques(mapped_adjacency_list, 0)[coordination_number[0]:coordination_number[1]+1]
    cliques = [np.vectorize(map2resind.get)(clique) for clique in cliques]

    #get neighbors and build flattened distance matrices
    for sub_cliques in cliques:
        for _clique in sub_cliques:
            clique = set(sum([_get_neighbors(structure, resind, no_neighbors) for resind in list(_clique)], []))
            selstr = 'resindex ' + ' '.join([str(i) for i in clique])
            clique_backbone = structure.select('protein').select('name N C CA O').select(selstr)

            padding = max_atoms - clique_backbone.numAtoms()
            flattened_dist_mat = np.lib.pad(buildDistMatrix(clique_backbone, clique_backbone), ((0,padding), (0,padding)), 'constant', constant_values=0).flatten()
            assert len(flattened_dist_mat) == 2304

            features.append(flattened_dist_mat)
            identifiers.append([resind2id[resind] for resind in clique])
            sources.append(pdb_file)

    site_df = pd.DataFrame({'features': features, 'identifiers': identifiers, 'sources': sources})
    return site_df