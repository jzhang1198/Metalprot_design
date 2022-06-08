#!/usr/bin/env python3

"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This script runs distributed prediction of metal binding.
"""

#imports 
import sys
import pickle
import numpy as np
from Metalprot_design.predictor import predict

def distribute_tasks(site_df: str):
    path2output = sys.argv[1] #path to store outputs  
    no_jobs = 1
    job_id = 0

    if len(sys.argv) > 3:
        no_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    row_indices = np.linspace(0, len(site_df)-1, len(site_df))
    task_rows = np.array_split(row_indices, no_jobs)[job_id]
    start_ind = int(task_rows[0])
    end_ind = int(task_rows[-1]) + 1
    tasks = site_df[start_ind:end_ind]
    print(f'Predicting coordinates for indices {start_ind}:{end_ind}')

    return path2output, job_id, tasks   

if __name__ == '__main__':
    SITE_DF_FILE = '/wynton/home/rotation/jzhang1198/data/metalprot_design/kehan_ion_channels/backbones/sites/compiled_site_dfs.pkl'

    with open(SITE_DF_FILE, 'rb') as f:
        site_df = pickle.load(f)

    path2output, job_id, tasks  = distribute_tasks(site_df)
    predict(path2output, job_id, tasks, path2model='./data/model_no_encodings')

