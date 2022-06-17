"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This file contains functions for loading structures for metal binding predictions.
"""

#imports 
from pypivoter.degeneracy_cliques import enumerateCliques
from prody import parsePDB, buildDistMatrix, AtomGroup
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import numpy as np
import pandas as pd
import itertools
import operator

def _get_helices(pdb_file: str, c_alphas: AtomGroup):
    """Identifies all helical residues and determines which helix they lie on.

    Args:
        pdb_file (str): Path to input pdb file.
        c_alphas (AtomGroup): An selection of all alpha carbons from the input pdb file.

    Returns:
        dict (dict): A dictionary that maps all helical residue indices to their corresponding helix index.
    """

    #compute dssp and find helical residues
    dssp_dict, _ = dssp_dict_from_pdb_file(pdb_file)
    helical_residues = [(resnum, chid) for resnum, chid in zip(c_alphas.getResnums(), c_alphas.getChids()) if dssp_dict[(chid, (' ', int(resnum), ' '))][1] in ['H', 'G', 'I']]
    helical_residues_by_chid = [list(group) for key, group in itertools.groupby(helical_residues, operator.itemgetter(1))]

    #identify residues that lie on a contiguous helix
    contiguous_helices = []
    for resis in helical_residues_by_chid:
        chid = resis[0][1]
        resnums = np.array([tup[0] for tup in resis])
        discontinuities = np.concatenate((np.array([0]), np.add(np.argwhere(resnums[1:] - resnums[0:-1] != 1).flatten(), 1), np.array([resnums[-1]])))
        contiguous_helices += [[(num, chid) for num in resnums[i:j]] for i,j in zip(discontinuities[0:-1], discontinuities[1:])]
    
    #create a dictionary mapping resindex to helix index
    dict = {}
    for ind in range(0, len(contiguous_helices)):
        for identifier in contiguous_helices[ind]:
            dict[c_alphas.select(f'chid {identifier[1]}').select(f'resnum {identifier[0]}').getResindices()[0]] = ind

    return dict

def _get_neighbors(structure, resind: int, no_neighbors: int):
    """Helper function. Finds neighbors of an input coordinating residue.

    Args:
        structure (prody.AtomGroup): Structure of full-length protein.
        coordinating_resind (int): Resindex of reference residue.
        no_neighbors (int): Number of neighbors in primary sequence to coordinating residue.

    Returns:
        core_fragment (list): List containing resindices of coordinating residue and neighbors. 
    """

    #find all resindices in the chain of input resind
    chain_id = list(set(structure.select(f'resindex {resind}').getChids()))[0]
    all_resinds = structure.select(f'chain {chain_id}').select('protein').getResindices()
    terminal = max(all_resinds)
    start = min(all_resinds) 

    #get fragment 
    extend = np.array(range(-no_neighbors, no_neighbors+1))
    _fragment = np.full((1,len(extend)), resind) + extend

    #if the input resind is a terminal residue, ensure that the non-existent resindex is removed
    fragment = [ind for ind in list(_fragment[ (_fragment >= start) & (_fragment <= terminal) ]) if ind in all_resinds] 

    return fragment

def _compute_ca_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray):
    """Computes imaginary alpha carbon - beta carbon bond vectors given N, CA, and C coordinates for all residues in a vectorized fashion.

    Args:
        n (np.ndarray): nx3 array containing N atom coordinates for all residues.
        ca (np.ndarray): nx3 array containing CA atom coordinates for all residues.
        c (np.ndarray): nx3 array containing C atom coordinates for all residues.

    Returns:
        ca_cb: nx3 array containing imaginary alpha carbon - beta carbon bond vectors.
    """

    #compute ca-n and ca-c bond vectors
    ca_n = np.vstack([(1 / np.linalg.norm(n - ca, axis=1))]*3).T * (n - ca)
    ca_c = np.vstack([(1 / np.linalg.norm(c - ca, axis=1))]*3).T * (c - ca)

    #using trigonometry, we can compute an imaginary ca-cb vector
    n1 = (ca_n + ca_c) * -1
    n2 = np.cross(ca_n, ca_c, axis=1)
    n1 = np.vstack([(1 / np.linalg.norm(n1, axis=1))]*3).T * n1
    n2 = np.vstack([(1 / np.linalg.norm(n2, axis=1))]*3).T * n2
    d = (1.54*np.sin(np.deg2rad(54.75))) * n2
    v = (1.54*np.cos(np.deg2rad(54.75))) * n1
    ca_cb = d+v

    return ca_cb

def _compute_angles(vec1: np.ndarray, vec2: np.ndarray):
    """Helper function for computing angles between bond vectors in a vectorized fashion.

    Args:
        vec1 (np.ndarray): An nx3 array containing bond vectors.
        vec2 (np.ndarray): Another nx3 array containing bond vectors.

    Returns:
        angles (np.ndarray): The angles between the vectors. 
    """

    dot = np.sum(vec1 * vec2, axis=1)
    norm1, norm2 = np.linalg.norm(vec1, axis=1), np.linalg.norm(vec2, axis=1)
    angles = np.degrees(np.arccos(dot / (norm1 * norm2)))
    return angles

def _filter_by_angle(edge_list: np.ndarray, structure: AtomGroup, distances: np.ndarray):
    """Filters pairs of contacts based on relative orientation of Ca-Cb and Ca-Ca bond vectors. 

    Args:
        edge_list (np.ndarray): nx2 array containing pairs of contacts.
        structure (AtomGroup): AtomGroup object of input structure.
        distances (np.ndarray): Array of length n containing distance between each contact.

    Returns:
        filtered (np.ndarray): nx2 array containing filtered contacts.
    """

    #get backbone atom coordinates for all residues included in the edge list
    all_resindices = set(np.concatenate(list(edge_list)))
    coordinates = dict([(resindex, structure.select('protein').select('name C CA N').select(f'resindex {resindex}').getCoords()) for resindex in all_resindices])

    #for each pair of contacts, get coordinates for atom i and j
    n_i, n_j = np.vstack([coordinates[resindex][0].flatten() for resindex in edge_list[:,0]]), np.vstack([coordinates[resindex][0].flatten() for resindex in edge_list[:,1]])
    ca_i, ca_j = np.vstack([coordinates[resindex][1].flatten() for resindex in edge_list[:,0]]), np.vstack([coordinates[resindex][1].flatten() for resindex in edge_list[:,1]])
    c_i, c_j = np.vstack([coordinates[resindex][2].flatten() for resindex in edge_list[:,0]]), np.vstack([coordinates[resindex][2].flatten() for resindex in edge_list[:,1]])

    #compute ca-cb bond vector for atom i and j
    ca_cb_i, ca_cb_j = _compute_ca_cb(n_i, ca_i, c_i), _compute_ca_cb(n_j, ca_j, c_j)

    #compute the ca-ca bond vector between atom i and j. also compute the ca-ca/ca-cbi and ca-ca/ca-bj angles.
    ca_i_ca_j = ca_j - ca_i
    angles_i, angles_j = _compute_angles(ca_cb_i, ca_i_ca_j), _compute_angles(ca_cb_j, ca_i_ca_j)

    #filter based on angle cutoffs
    accepted = np.argwhere(distances <= 7)
    filtered_inds = np.intersect1d(np.intersect1d(np.argwhere(angles_i < 130), np.argwhere(angles_j > 30)), np.argwhere(distances > 7))
    filtered = edge_list[np.union1d(accepted, filtered_inds)]
    return filtered

def _filter_by_no_helices(cliques: list, mappings: dict, helix_cutoff: int):
    """Filters cliques based on number of participating helices.

    Args:
        cliques (list): Cliques computed by enumerateCliques.
        mappings (dict): Dictionary that maps resindices to helix indices.
        helix_cutoff (int): Minimum mumber of helices within a putative binding core.

    Returns:
        filtered_cliques (list): Cliques satisfying helix cutoff.
    """

    filtered_cliques = []
    for clique in cliques:
        #map resindices to helix index in a vectorized fashion
        mapped = np.vectorize(mappings.get)(clique)
        mapped = np.sort(mapped, axis=1)

        #count number of helices and get indices of cliques that satisfy cutoff
        no_helices = (mapped[:,1:] != mapped[:,:-1]).sum(axis=1) + 1
        indices = np.argwhere(no_helices >= helix_cutoff)
        filtered_cliques.append(clique[indices].squeeze())

    return filtered_cliques

def identify_sites_rational(pdb_file: str, cutoff: float, helix_cutoff: int, coordination_number=(2,4), no_neighbors=1):
    """Main function that runs site enumeration.

    Args:
        pdb_file (str): Path to input pdb file.
        cutoff (float): Defines upper limit on Ca distance between binding core residues.
        helix_cutoff (int): Defines lower limit on number of helices in a putative binding core.
        coordination_number (tuple, optional): Defines lower and upper limit on coordination numbers considered. Defaults to (2,4).
        no_neighbors (int, optional): Defines number of neighbors in primary sequence to include for clique residues. Defaults to 1.

    Returns:
        site_df (pd.DataFrame): Contains flattened distance matrices, residue numbers and chain IDs, and source pdb file paths for all enumerated cores.
    """

    features, identifiers, sources = [], [], []
    max_atoms = 4 * (coordination_number[1] + (coordination_number[1] * no_neighbors * 2))
    structure = parsePDB(pdb_file)
    c_alphas = structure.select('protein').select('name CA')

    #enumerate all helical residues in the input structure
    helix_map = _get_helices(pdb_file, c_alphas)
    helix_resindices = list(helix_map.keys())
    helix_resindices.sort()
    helix_resindices_selstr = 'resindex ' + ' '.join([str(i) for i in helix_resindices])

    #create dictionary to map resindices to chain IDs and residue numbers
    resind2id = dict([(resindex, (structure.select(f'resindex {resindex}').getResnums()[0], structure.select(f'resindex {resindex}').getChids()[0])) for resindex in set(structure.select('protein').getResindices())])

    #build alpha carbon distance matrix. iterate through columns and enumerate pairs of residues that are within a cutoff distance.
    ca_dist_mat = buildDistMatrix(c_alphas.select(helix_resindices_selstr), c_alphas.select(helix_resindices_selstr))
    edge_list = []
    edge_weights = np.array([])
    row_indexer = 0
    for col_ind in range(len(ca_dist_mat)):
        for row_ind in range(1+row_indexer, len(ca_dist_mat)):
            distance = ca_dist_mat[row_ind, col_ind]

            if distance <= cutoff:
                edge_list.append(np.array([helix_resindices[col_ind], helix_resindices[row_ind]]))
                edge_weights = np.append(edge_weights, distance)

        row_indexer += 1

    #filter by relative orientation of ca-cb bond vectors
    edge_list = _filter_by_angle(np.vstack(edge_list), structure, edge_weights)
    resind2map = dict([(resind, index) for index, resind in enumerate(np.sort(np.unique(edge_list.flatten())))])
    map2resind = dict([(index, resind) for index, resind in enumerate(np.sort(np.unique(edge_list.flatten())))])
    mapped_edge_list = np.vectorize(resind2map.get)(edge_list)
    assert np.max(mapped_edge_list) == len(np.unique(mapped_edge_list)) - 1

    #enumerate cliques
    cliques = enumerateCliques(mapped_edge_list, coordination_number[1])[coordination_number[0]:]
    cliques = [np.vectorize(map2resind.get)(clique) for clique in cliques]

    #if a helix cutoff is provided, filter by number of donating helices
    if helix_cutoff:
        cliques = _filter_by_no_helices(cliques, helix_map, helix_cutoff)

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