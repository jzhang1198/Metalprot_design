"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This file contains functions for predicting metal binding properties of a given protein.
"""

#imports
import os
import json
import torch
import numpy as np
import pandas as pd
from prody import parsePDB
from Metalprot_design.model import DoubleLayerNet

def _load_neural_net(path2model: str):
    weights_file = os.path.join(path2model, 'model.pth')
    with open(os.path.join(path2model, 'config.json'), 'r') as f:
        config = json.load(f)

    model = DoubleLayerNet(config['input'], config['l1'], config['l2'], config['l3'], config['output'], config['input_dropout'], config['hidden_dropout'])
    model.load_state_dict(torch.load(weights_file, map_location='cpu'))
    model.eval()
    return model

def _triangulate(backbone_coords, distance_prediction):
    distance_prediction = distance_prediction[0:len(backbone_coords)]

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
    
    result = scipy.optimize.minimize(objective, guess)
    solution = result.x
    rmsd = objective(solution)

    return solution, rmsd

def _extract_coordinates(source_file: str, identifier_permutation):
    """_summary_

    Args:
        source_file (str): _description_
        positive (bool, optional): _description_. Defaults to False.
    """

    core = parsePDB(source_file)
    for iteration, id in enumerate(identifier_permutation):
        residue = core.select(f'chain {id[1]}').select(f'resnum {id[0]}').select('name C CA N O').getCoords()
        coordinates = residue if iteration == 0 else np.vstack([coordinates, residue])

    return coordinates

def predict(path2output: str, job_id: int, site_df: pd.DataFrame, path2model: str):

    model = _load_neural_net(path2model)
    X = np.vstack([array for array in site_df['features']])
    prediction = model.forward(torch.from_numpy(X)).cpu().detach().numpy()

    deviation = np.array([np.nan] * len(prediction))
    completed = 0
    for distance_prediction, pointer, resindex_permutation in zip(prediction, list(site_df['source']), list(site_df['identifiers'])):
        try:
            source_coordinates = _extract_coordinates(pointer, resindex_permutation)
            solution, rmsd = _triangulate(source_coordinates, distance_prediction)
            completed += 1

        except:
            solution, rmsd = np.array([np.nan, np.nan, np.nan]), np.nan

        if 'solutions' not in locals():
            solutions = solution
            rmsds = rmsd

        else:
            solutions = np.vstack([solutions, solution])
            rmsds = np.append(rmsds, rmsd)

    predictions = pd.DataFrame({'predicted_distances': list(prediction),
        'predicted_coordinates': list(solutions),
        'confidence': rmsds,
        'deviation': deviation,
        'barcode': site_df['barcode'].to_numpy()})

    predictions.to_pickle(os.path.join(path2output, f'predictions{job_id}.pkl'))

    print(f'Coordinates and RMSDs computed for {completed} out of {len(prediction)} observations.')