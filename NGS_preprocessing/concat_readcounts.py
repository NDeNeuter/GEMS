#!/usr/bin/env python3

import os

import pandas as pd

def read_df(file):
    
    df = pd.read_csv(file, sep='\t')
    df.set_index('genename', inplace=True)
    df.fillna(0, inplace=True)
    return df


raw_readcounts = []
for dr, subdr, files in os.walk('..'):
    for file in files:
        if file.endswith('readcounts.txt'):
            raw_readcounts.append('{}/{}'.format(dr, file))
            
to_concat = []
for file in raw_readcounts:
    df = read_df(file)
    # extract subject name from column names to group on
    df = df.T
    df['name'] = df.index
    df['sample'] = df['name'].apply(lambda x: '_'.join(x.split('/')[-1].split('_')[:-3]))
    # group on subject name: puts all reads from same subject (but different lanes) together
    grouped_df = df.groupby('sample').sum().T
    # add df to list of dfs to concatenate together
    to_concat.append(grouped_df)
    
# concat dataframes together
concat_df = pd.concat(to_concat, axis=1).fillna(0)
# write df to file
concat_df.to_csv('../data/combined_readcounts.tsv', sep='\t')
