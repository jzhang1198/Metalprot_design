#!/usr/bin/env python3

"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This script runs distributed site identification on an input set of pdb files.
"""

#imports
import os
import sys
import pandas as pd
from Metalprot_design.loader import identify_sites_rational

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

    PDB_DIR = '/wynton/home/rotation/jzhang1198/scratch/test_enumeartion' #directory containing pdb files
    CUTOFF = 12 #cutoff distance (in angstroms) that defines a site. i've found that 12A should do the trick for almost all metal binding sites

    tasks, path2output, job_id = distribute_tasks(PDB_DIR)
    for task in tasks:
        name = task.split('/')[-1].split('.')[0]
        site_df = identify_sites_rational(task, CUTOFF) #identify putative binding cores and get distance matrices
        site_df.to_pickle(os.path.join(path2output, f'{name}_site_df{job_id}.pkl')) #write dataframe to a pickle file