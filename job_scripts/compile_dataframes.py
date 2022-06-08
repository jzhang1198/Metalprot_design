#!/usr/bin/env python3

"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This script stacks dataframes outputted by a parallel run into a single dataframe.
"""

#imports
import os
import sys
import pickle
import pandas as pd

def remove_duplicates(df):
    identifiers = df['identifiers']
    df.drop('identifiers', inplace=True, axis=1)
    identifiers = [tuple(x) for x in list(identifiers)]
    df['identifiers'] = identifiers
    df.drop_duplicates(subset=['identifiers', 'sources'], keep='first', inplace=True)
    
    return df

if __name__ == '__main__':
    path2output = sys.argv[1]
    df_files = [os.path.join(path2output, file) for file in os.listdir(path2output) if '.pkl' in file]

    total = len(df_files)
    counter = 0
    for file in df_files:
        with open(file, 'rb') as f:
            _df = pickle.load(f)

        df = pd.DataFrame(columns=_df.columns)
        df = pd.concat([df, _df])

        if len(df) > 222500:
            df = remove_duplicates(df)
            df.to_pickle(os.path.join(path2output, f'compiled_site_dfs{counter}.pkl'))
            df = pd.DataFrame(columns=_df.columns)

    df = remove_duplicates(df)
    df.to_pickle(os.path.join(path2output, f'compiled_site_dfs{counter}.pkl'))