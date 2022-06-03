"""
Author: Jonathan Zhang <jon.zhang@ucsf.edu>

This script stacks dataframes outputted by a parallel run into a single dataframe.
"""

#imports
import os
import sys
import pickle
import pandas as pd

if __name__ == '__main__':
    path2output = sys.argv[1]
    df_files = [os.path.join(path2output, file) for file in os.listdir(path2output) if 'site_df' and '.pkl' in file]

    for file in df_files:
        with open(file, 'rb') as f:
            _df = pickle.load(f)

        df = pd.concat([df, _df]) if 'df' in locals() else _df

    df.to_pickle(os.path.join(path2output, 'compiled_site_dfs.pkl'))