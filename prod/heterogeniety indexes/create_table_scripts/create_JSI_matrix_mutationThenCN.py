# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:13:03 2020

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
unique_mutation_columns = ['Patient ID','REF','ALT','GENE','EFFECT','CHROM','POSITION']

# =============================================================================
# Helpers
# =============================================================================

def get_mutation_jsi(muts, samples):
    muts_tmp = muts[muts['Sample ID'].isin(samples)].copy()
    muts_unique = muts_tmp.drop_duplicates(unique_mutation_columns)[unique_mutation_columns]
    muts_unique = muts_unique.set_index(unique_mutation_columns)
    muts_unique['Sample1'] = 0
    muts_unique['Sample2'] = 0
    muts_tmp = muts_tmp.set_index(unique_mutation_columns)
    muts_sample1 = muts_tmp[muts_tmp['Sample ID'] == samples[0]]
    muts_sample2 = muts_tmp[muts_tmp['Sample ID'] == samples[1]]
    for index, row in muts_unique.iterrows():
        if index in muts_sample1.index:
            muts_unique.at[index, 'Sample1'] = 1
        if index in muts_sample2.index:
            muts_unique.at[index, 'Sample2'] = 1
    muts_unique['Shared'] = 1
    muts_unique.loc[muts_unique['Sample1'] != muts_unique['Sample2'], 'Shared'] = 0
    return muts_unique[['Shared']]

def get_cn_jsi(cn,samples):
    cn_tmp = cn[cn['Sample ID'].isin(samples)].copy()
    cn_tmp['Copy_num'] = cn_tmp['Copy_num'].replace({2:1,-2:-1})
    cn_tmp = cn_tmp.pivot(index = 'GENE', columns = 'Sample ID', values = 'Copy_num')
    cn_tmp = cn_tmp[(cn_tmp[samples[0]] != 0)|(cn_tmp[samples[1]] != 0)]
    cn_tmp['Shared'] = 1
    cn_tmp.loc[cn_tmp[samples[0]] != cn_tmp[samples[1]], 'Shared'] = 0
    return cn_tmp[['Shared']]

# =============================================================================
# Import mutations, CN, TC, and clinical data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/%s_mutations_inclDependent.tsv' % cohort, sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/%s_cna.tsv' % cohort, sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc = tc[tc['Cohort'] == cohort]
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100
tc = tc[tc['Final tNGS_TC'] > 0]

tc_dict = dict(zip(tc['Sample ID'].tolist(), tc['Final tNGS_TC'].tolist()))

# =============================================================================
# Create same patient sample combination list
# =============================================================================

ls = []

for pt in tc['Patient ID'].unique().tolist():
    samples = tc[tc['Patient ID'] == pt]['Sample ID'].tolist()
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
matrix = pd.DataFrame(columns = tc['Sample ID'].unique(),index = tc['Sample ID'].unique())

for index, row in jsi.iterrows():
    sample1 = index[0]
    sample2 = index[1]
    mutation_jsi = get_mutation_jsi(muts, [sample1,sample2])
    cn_jsi = get_cn_jsi(cn,[sample1,sample2])
    if len(mutation_jsi) > 0:
        if mutation_jsi['Shared'].sum() == len(mutation_jsi):
            if len(cn_jsi) > 0:
                value = cn_jsi['Shared'].sum() / len(cn_jsi)
            else: 
                value = 1
        else:
            value = mutation_jsi['Shared'].sum() / len(mutation_jsi)
    elif len(cn_jsi) > 0:
        value = cn_jsi['Shared'].sum() / len(cn_jsi)
    else:
        value = 0

    matrix.at[sample1, sample2] = value
    print(value)


matrix.to_csv('G:\Andy Murtha\Ghent\M1RP\prod\heterogeniety indexes\jsi_matrix_mutThenCN.tsv', sep = '\t')