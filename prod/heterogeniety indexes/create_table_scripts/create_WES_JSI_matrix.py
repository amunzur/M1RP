# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 12:26:58 2021

@author: amurtha
"""

import pandas as pd
import numpy as np
import itertools as it
import multiprocessing as mp
from joblib import Parallel, delayed

# =============================================================================
# Import mutations
# =============================================================================

# called = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/june2021_wes_mutations_betastasis_melted_unfiltered.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

if 'muts' not in locals():
    muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/betastasis/june2021_wes_mutations_betastasis.tsv', sep = '\t', index_col=None)

# =============================================================================
# Clean data
# =============================================================================

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc.set_index('Sample ID')

muts = muts[[col for col in muts.columns.tolist() if 'WBC' not in col and 'gDNA' not in col and 'NL' not in col]]

# =============================================================================
# 
# =============================================================================

samples = muts.columns.tolist()[7:]
index_cols = ['CHROM', 'POSITION', 'GENE', 'REF','ALT', 'EFFECT']


def get_jsi(s1, s2, muts):
    if s1 == s2:
        return (1,None);
    s_muts = muts[index_cols+[s1,s2]].copy()
    s_muts = s_muts[(s_muts[s1].str.contains("*", regex = False))|(s_muts[s2].str.contains("*", regex = False))]
    # Get the mutant allele frequency of each mutation
    s_muts[s1+'_af'] = s_muts[s1].str.split('%').str[0].astype(float) / 100
    s_muts[s2+'_af'] = s_muts[s2].str.split('%').str[0].astype(float) / 100
    # Get depth at each mutation
    s_muts[s1+'_depth'] = s_muts[s1].str.split('(').str[1].str.split(')').str[0].astype(float)
    s_muts[s2+'_depth'] = s_muts[s2].str.split('(').str[1].str.split(')').str[0].astype(float)
    # Get the number of mutant reads for each mutation
    s_muts[s1+'_mutant_reads'] = s_muts[s1+'_af']*s_muts[s1+'_depth']
    s_muts[s2+'_mutant_reads'] = s_muts[s2+'_af']*s_muts[s2+'_depth']
    # Add tumor content as column
    s_muts[s1+'_tc'] = tc.at[s1,'Final tNGS_TC']
    s_muts[s2+'_tc'] = tc.at[s2,'Final tNGS_TC']
    # Calculate expected number of mutant reads per mutation
    s_muts[s1+'_expected_mReads'] = s_muts[s1+'_tc'] * s_muts[s1+'_depth'] / 2
    s_muts[s2+'_expected_mReads'] = s_muts[s2+'_tc'] * s_muts[s2+'_depth'] / 2
    # Add columns for called in each sample
    s_muts[s1+'_called'] = s_muts[s1].str.contains('*', regex = False)
    s_muts[s2+'_called'] = s_muts[s2].str.contains('*', regex = False)
    # set called to true if MAF > 0 and mutant reads > 3
    s_muts.loc[(s_muts[s1+'_af'] > 0) & (s_muts[s1+'_mutant_reads'] >= 3), s1+'_called'] = True
    s_muts.loc[(s_muts[s2+'_af'] > 0) & (s_muts[s2+'_mutant_reads'] >= 3), s2+'_called'] = True
    # Remove rows where expected mutant reads is < 3 in either sample
    s_muts = s_muts[(s_muts[s1+'_expected_mReads'] >= 3)&(s_muts[s2+'_expected_mReads'] >= 3)]
    # Return in no mutations are present
    if len(s_muts) == 0:
        return (np.nan, 0)
    # Add shared column. Calculate JSI from this column and return value
    s_muts['shared'] = (s_muts[s1+'_called'] == True)&(s_muts[s2+'_called'] == True)
    s_muts['shared'] = s_muts['shared'].map({True:1,False:0})
    return (s_muts['shared'].sum() / len(s_muts), len(s_muts))

matrix = pd.DataFrame(index = samples, columns = samples)
matrix_n = pd.DataFrame(index = samples, columns = samples)

def func(s1, s2):
    return get_jsi(s1, s2, muts)

results = [func(s1,s2) for (s1,s2) in it.combinations_with_replacement(samples, 2)]

n_cores = mp.cpu_count() 
results = Parallel(n_jobs=n_cores)(delayed(func)(s1,s2) for (s1,s2) in it.combinations_with_replacement(samples, 2))

for (s1,s2),r in zip(it.combinations_with_replacement(samples, 2), results):
    matrix.at[s1,s2] = r[0]
    matrix.at[s2,s1] = r[0]
    matrix_n.at[s1,s2] = r[1]
    matrix_n.at[s2,s1] = r[1]
    
matrix.to_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES.tsv', sep = '\t')
matrix_n.to_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_nMuts.tsv', sep = '\t')
    