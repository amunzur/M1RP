# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:41:58 2020

@author: amurtha
"""

import seaborn as sns
import matplotlib.gridspec
import pandas as pd
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Constants
# =============================================================================

cohort = 'M1RP'
unique_mutation_columns = ['Patient ID','REF','ALT','GENE','EFFECT','CHROM',]

# =============================================================================
# Helpers
# =============================================================================

def get_mutation_jsi(muts, samples):
    muts = muts.sort_index()
    muts_tmp = muts[muts['Sample ID'].isin(samples)].copy()
    muts_unique = muts_tmp.drop_duplicates(unique_mutation_columns)[unique_mutation_columns]
    muts_unique = muts_unique.set_index(unique_mutation_columns)
    muts_unique['Sample1'] = 0
    muts_unique['Sample2'] = 0

    muts_tmp = muts_tmp.set_index(unique_mutation_columns)
    muts_sample1 = muts_tmp[muts_tmp['Sample ID'] == samples[0]]
    muts_sample2 = muts_tmp[muts_tmp['Sample ID'] == samples[1]]
    muts_sample1 = muts_sample1.sort_index()
    muts_sample2 = muts_sample2.sort_index()
    muts_unique = muts_unique.sort_index()
    for index, row in muts_unique.iterrows():
        if index in muts_sample1.index:
            muts_unique.at[index, 'Sample1'] = 1
        if index in muts_sample2.index:
            muts_unique.at[index, 'Sample2'] = 1
    muts_unique['Shared'] = 1
    muts_unique.loc[muts_unique['Sample1'] != muts_unique['Sample2'], 'Shared'] = 0
    return muts_unique[['Shared']]

# =============================================================================
# Import mutations, CN, TC, and clinical data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/wxs_muts_all.tsv', sep = '\t')

# =============================================================================
# Create same patient sample combination list
# =============================================================================

ls = []

for pt in muts['Patient ID'].unique().tolist():
    samples = muts[muts['Patient ID'] == pt]['Sample ID'].unique().tolist()
    ls.append(list(itertools.product(samples, samples)))
    
samples = []
for l in ls:
    samples = samples+l
    
del pt, ls

samples = [sample for sample in samples if sample[0] != sample[1]]
index = pd.MultiIndex.from_tuples(samples, names=['Sample1', 'Sample2'])

# =============================================================================
# Create dataframe
# =============================================================================

jsi = pd.DataFrame(index = index)
matrix = pd.DataFrame(columns = muts['Sample ID'].unique(),index = muts['Sample ID'].unique())

for index, row in jsi.iterrows():
    sample1 = index[0]
    sample2 = index[1]
    mutation_jsi = get_mutation_jsi(muts, [sample1,sample2])
    # cn_jsi = get_cn_jsi(cn,[sample1,sample2])
    if len(mutation_jsi) > 0:
        value = mutation_jsi['Shared'].sum() / len(mutation_jsi)
    else:
        value = 1
    matrix.at[sample1, sample2] = value


matrix.to_csv('G:\Andy Murtha\Ghent\M1RP\dev\heterogeniety indexes\jsi_matrix_wxs.tsv', sep = '\t')