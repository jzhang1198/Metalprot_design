"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This file contains functions for predicting metal binding properties of a given protein.
"""

#imports
import os
import json
import torch
import scipy
import numpy as np
import pandas as pd
from prody import parsePDB, AtomGroup
from Metalprot_design.model import DoubleLayerNet

def _load_neural_net(path2model: str):
    """Helper function that loads the trained regressor.

    Args:
        path2model (str): The path to the directory containing the model weights and config file.

    Returns:
        model: Loaded model.
    """
    weights_file = os.path.join(path2model, 'model.pth')
    with open(os.path.join(path2model, 'config.json'), 'r') as f:
        config = json.load(f)

    model = DoubleLayerNet(config['input'], config['l1'], config['l2'], config['l3'], config['output'], config['input_dropout'], config['hidden_dropout'])
    model.load_state_dict(torch.load(weights_file, map_location='cpu'))
    model.eval()
    return model

def _extract_coordinates(core: AtomGroup, identifiers: list):
    """Helper function for loading coordinates

    Args:a
        core (AtomGroup): Atom group object of a source pdb structure.
        identifiers (list): List of tuples, where the elements of the tuple are a residue number and a chain ID corresponding to a core residue.

    Returns:
        coordinates (np.ndarray): nx3 array containing backbone coordinates of the core residues.
    """
    
    resinds = [core.select(f'chid {tup[1]}').select(f'resnum {tup[0]}').getResindices()[0] for tup in identifiers]
    selstr = 'resindex ' + ' '.join([str(i) for i in resinds])
    coordinates = core.select(selstr).getCoords() 

    return coordinates

def _triangulate(backbone_coords: np.ndarray, distance_prediction: np.ndarray):
    """Given coordinates and distances, computes the metal coordinates via triangulation.

    Args:
        backbone_coords (np.ndarray): nx3 array containing backbone atom coordinates.
        distance_prediction (np.ndarray): Array of length n containing backbone atom-metal distances.

    Returns:
        solution (np.ndarray): Array of length 3 containing x, y, and z coordinates of the metal.
        rmsd (float): RMSD between the predicted distances and the distances realized by placement of the metal in coordinate space.
    """

    distance_prediction = distance_prediction[0:len(backbone_coords)]

    #define objective function and initial guess
    guess = backbone_coords[0]
    def objective(v):
        x,y,z = v
        distances = np.zeros(backbone_coords.shape[0])
        for i in range(0, backbone_coords.shape[0]):
            atom = backbone_coords[i]
            dist = np.linalg.norm(atom - np.array([x,y,z]))
            distances[i] = dist
        rmsd = np.sqrt(np.mean(np.square(distances - distance_prediction)))
        return rmsd
    
    #compute coordinates via minimization of objective
    result = scipy.optimize.minimize(objective, guess)
    solution = result.x
    rmsd = objective(solution)

    return solution, rmsd

def predict(path2output: str, job_id: int, site_df: pd.DataFrame, path2model: str):
    """Main function for prediction of metal coordinates.

    Args:
        path2output (str): Path to output directory to write files to.
        job_id (int): The id for a given task.
        site_df (pd.DataFrame): Dataframe containing distance matrices, sources, and identifiers for enumerated cores.
        path2model (str): Path to the directory containing regressor weights and config file.
    """

    #load model and compute a forward pass through the regressor
    model = _load_neural_net(path2model)
    X = np.vstack([array for array in site_df['features']])
    prediction = model.forward(torch.from_numpy(X)).cpu().detach().numpy()

    #generate a dictionary that maps coordinates to the corresponding source file
    ag_dict = dict([(pointer, parsePDB(pointer).select('protein').select('name N C CA O')) for pointer in set(list(site_df['sources']))])

    completed = 0
    for distance_prediction, pointer, identifiers in zip(prediction, list(site_df['sources']), list(site_df['identifiers'])):

        #compute metal coordinates and rmsds
        try:
            source_coordinates = _extract_coordinates(ag_dict[pointer], identifiers)
            solution, rmsd = _triangulate(source_coordinates, distance_prediction)
            completed += 1

        except:
            solution, rmsd = np.array([np.nan, np.nan, np.nan]), np.nan

        #append coordinates and rmsds to a vector for tracking
        if 'solutions' not in locals():
            solutions = solution
            rmsds = rmsd

        else:
            solutions = np.vstack([solutions, solution])
            rmsds = np.append(rmsds, rmsd)

    #write output prediction file
    predictions = pd.DataFrame({'predicted_distances': list(prediction),
        'predicted_coordinates': list(solutions),
        'confidence': rmsds,
        'sources': list(site_df['sources']),
        'identifiers': list(site_df['identifiers'])})
    predictions.to_pickle(os.path.join(path2output, f'predictions{job_id}.pkl'))

    print(f'Coordinates and RMSDs computed for {completed} out of {len(prediction)} observations.')