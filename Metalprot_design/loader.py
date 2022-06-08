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
    dssp_dict, _ = dssp_dict_from_pdb_file(pdb_file)
    helical_residues = [(resnum, chid) for resnum, chid in zip(c_alphas.getResnums(), c_alphas.getChids()) if dssp_dict[(chid, (' ', int(resnum), ' '))][1] in ['H', 'G', 'I']]
    helical_residues_by_chid = [list(group) for key, group in itertools.groupby(helical_residues,operator.itemgetter(1))]

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

def _identify_subclusters(identifiers: list, contiguous_helices: dict, coordination_number: tuple):
    subclusters = []
    for number in range(coordination_number[0], coordination_number[1]+1):
        all_combinations = list(itertools.combinations(identifiers, number))
        no_non_degenerate_helices = np.array([len(set(contiguous_helices[tup] for tup in combination)) for combination in all_combinations])
        candidate_combinations = [all_combinations[int(ind)] for ind in np.argwhere(no_non_degenerate_helices > 1).flatten()]
        subclusters += candidate_combinations

    return subclusters

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

def identify_sites(pdb_file: str, cuttoff: float, coordination_number: tuple, no_neighbors: int):

    features, identifiers, sources, barcodes = [], [], [], []

    max_atoms = 4 * (coordination_number[1] + (coordination_number[1] * no_neighbors * 2))
    structure = parsePDB(pdb_file)
    c_alphas = structure.select('protein').select('name CA')

    contiguous_helices = _get_helices(pdb_file, c_alphas)
    helix_resindices = [structure.select(f'chid {tup[1]}').select(f'resnum {tup[0]}').getResindices()[0] for tup in contiguous_helices.keys()]
    helix_resindices_selstr = ' '.join([str(i) for i in helix_resindices])

    counter = 0
    for resnum, chid in  contiguous_helices.keys():
        resindex = c_alphas.select(f'chid {chid}').select(f'resnum {resnum}').getResindices()[0]
        c_alpha_cluster = c_alphas.select(f'within {cuttoff} of resindex {resindex}').select(f'resindex {helix_resindices_selstr}')

        for subcluster in _identify_subclusters([(num, chid) for num, chid in zip(c_alpha_cluster.getResnums(), c_alpha_cluster.getChids())], contiguous_helices, coordination_number):
            cluster_resindices = np.array([structure.select(f'protein').select(f'chid {tup[1]}').select(f'resnum {tup[0]}').getResindices()[0] for tup in subcluster])
            cluster_resindices = sum([_get_neighbors(structure, resind, no_neighbors) for resind in cluster_resindices], [])
            selstr = 'resindex ' + ' '.join([str(i) for i in cluster_resindices])
            core_backbone = structure.select('protein and name CA C N O').select(selstr)

            padding = max_atoms - core_backbone.numAtoms()
            flattened_dist_mat = np.lib.pad(buildDistMatrix(core_backbone, core_backbone), ((0,padding), (0,padding)), 'constant', constant_values=0).flatten()
            assert len(flattened_dist_mat) == 2304

            features.append(flattened_dist_mat)
            identifiers.append([(resnum, chain) for resnum, chain in zip(core_backbone.select('name CA').getResnums(), core_backbone.select('name CA').getChids())])
            sources.append(pdb_file)
            barcodes.append(counter)

            counter += 1

        else:
            continue

    site_df = pd.DataFrame({'features': features, 'identifiers': identifiers, 'sources': sources, 'barcodes': barcodes})
    return site_df


def _get_ca_cb(backbone_coords: np.ndarray):
    ca_n = backbone_coords[0] - backbone_coords[1]
    ca_n = ca_n / np.linalg.norm(ca_n)
    ca_c = backbone_coords[2] - backbone_coords[1]
    ca_c = ca_c / np.linalg.norm(ca_c)
    n2 = np.cross(ca_n, ca_c) 
    n2 = n2 / np.linalg.norm(n2)
    n1 = ((ca_n + ca_c) / np.linalg.norm(ca_n + ca_c)) * -1
    
    d = (1.54*np.sin(np.deg2rad(54.75)))*n2
    v = (1.54*np.cos(np.deg2rad(54.75)))*n1

    ca_cb = d + v
    ori = backbone_coords[1]

    return ca_cb, ori

def _filter(adjacency_list: np.ndarray, structure: AtomGroup, angle_cutoff: float):

    filtered = []
    coords_i = structure.select('protein').select(f'resindex {adjacency_list[0,0]}').select('name C CA N O').getCoords()
    for adjacency in adjacency_list:
        coords_j =  structure.select('protein').select(f'resindex {adjacency[1]}').select('name C CA N O').getCoords()
        cb_veci = _get_ca_cb(coords_i)
        cb_vecj = _get_ca_cb(coords_j)

        dot = cb_veci.dot(cb_vecj)
        angle = np.degrees(np.arccos(dot / (np.linalg.norm(cb_veci) * np.linalg.norm(cb_vecj))))
        
        if angle < angle_cutoff:
            filtered.append(adjacency)

    return filtered

def identify_sites_rational(pdb_file: str, cuttoff: float, angle_cutoff: float, coordination_number: tuple, no_neighbors: int):
    features, identifiers, sources = [], [], []

    max_atoms = 4 * (coordination_number[1] + (coordination_number[1] * no_neighbors * 2))
    structure = parsePDB(pdb_file)
    c_alphas = structure.select('protein').select('name CA')

    #enumerate all helical residues in the input structure
    contiguous_helices = _get_helices(pdb_file, c_alphas)
    helix_resindices = [structure.select(f'chid {tup[1]}').select(f'resnum {tup[0]}').getResindices()[0] for tup in contiguous_helices.keys()]
    helix_resindices.sort()
    helix_resindices_selstr = 'resindex ' + ' '.join([str(i) for i in helix_resindices])
    helix_identifiers = dict([(resindex, (structure.select(f'resindex {resindex}').getResnums()[0], structure.select(f'resindex {resindex}').getChids()[0])) for resindex in helix_resindices])

    #build alpha carbon distance matrix. iterate through columns and enumerate pairs of residues that are within a cutoff distance.
    ca_dist_mat = buildDistMatrix(c_alphas.select(helix_resindices_selstr), c_alphas.select(helix_resindices_selstr))
    adjacency_list = []
    row_indexer = 0
    for col_ind in range(len(ca_dist_mat)):
        for row_ind in range(1+row_indexer, len(ca_dist_mat)):
            distance = ca_dist_mat[row_ind, col_ind]

            if distance <= cuttoff:
                adjacency_list.append(np.array([helix_resindices[col_ind], helix_resindices[row_ind]]))

        #remove pairs of residues that have diverging alpha caarbon-beta carbon vectors
        filtered_adjacency_list = np.vstack(_filter(adjacency_list, structure, angle_cutoff))
        cliques = enumerateCliques(filtered_adjacency_list, coordination_number[1])
        
        #drop cliques that are under lower limit on coordination number

        #get neighbors and build flattened distance matrices
        for _clique in cliques:
            clique = sum([_get_neighbors(structure, resind, no_neighbors) for resind in list(_clique)], [])
            selstr = 'resindex ' + ' '.join([str(i) for i in clique])
            clique_backbone = structure.select('protein').select('name N C CA O').select(selstr)

            padding = max_atoms - clique_backbone.numAtoms()
            flattened_dist_mat = np.lib.pad(buildDistMatrix(clique_backbone, clique_backbone), ((0,padding), (0,padding)), 'constant', constant_values=0).flatten()
            assert len(flattened_dist_mat) == 2304

            features.append(flattened_dist_mat)
            identifiers.append([helix_identifiers[resind] for resind in clique])
            sources.append(pdb_file)

        row_indexer += 1

    site_df = pd.DataFrame({'features': features, 'identifiers': identifiers, 'sources': sources})
    return site_df