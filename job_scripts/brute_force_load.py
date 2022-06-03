#!/usr/bin/env python3

"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This script runs distributed site identification on an input set of pdb files.
"""

#imports
import os
import sys
import pandas as pd
from Metalprot_design.loader import identify_sites

def distribute_tasks(pdb_dir: str):
    path2output = sys.argv[1] #path to store outputs  
    no_jobs = 1
    job_id = 0

    if len(sys.argv) > 3:
        no_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    all_tasks = [os.path.join(pdb_dir, file) for file in os.listdir(pdb_dir) if 'pdb' in file]
    tasks = [all_tasks[i] for i in range(0, len(all_tasks)) if i % no_jobs == job_id]

    return tasks, path2output, job_id

if __name__ == '__main__':
    PDB_DIR = '/Users/jonathanzhang/Documents/scratch/test' #directory containing pdb files
    CUTOFF = 10 #cutoff distance (in angstroms) that defines a site
    COORDINATION_NUMBER = (2,4) #don't change bottom two variables. the model was trained with 2-4 coordinate binding sites including the neighrbors to the immediate left or right of a given coordinating residue
    NO_NEIGHBORS = 1

    tasks, path2output, job_id = distribute_tasks(PDB_DIR)
    for task in tasks:
        # try:
        _site_df = identify_sites(task, CUTOFF, COORDINATION_NUMBER, NO_NEIGHBORS) #identify putative binding cores and get distance matrices
        site_df = pd.concat([site_df, _site_df]) if 'site_df' in locals() else _site_df

        # except:
        #     with open(os.path.join(path2output, 'failed.txt'), 'a') as f:
        #         f.write(task + '\n')

    site_df.to_pickle(os.path.join(path2output, f'site_df{job_id}.pkl')) #write dataframe to a pickle file