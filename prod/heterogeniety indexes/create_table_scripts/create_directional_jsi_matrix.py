# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:53:09 2021

@author: amurtha
"""

import pandas as pd
import numpy as np
import itertools as it
import multiprocessing as mp
from joblib import Parallel, delayed
import string

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

samples = muts.columns.tolist()[8:]
index_cols = ['CHROM', 'POSITION', 'GENE', 'REF','ALT', 'EFFECT']


def get_p2m_jsi(sp, sm, muts):
    s_muts = muts[index_cols+[sp,sm]].copy()
    s_muts = s_muts[(s_muts[sp].str.contains("*", regex = False))]
    # Get the mutant allele frequency of each mutation
    s_muts[sp+'_af'] = s_muts[sp].str.split('%').str[0].astype(float) / 100
    s_muts[sm+'_af'] = s_muts[sm].str.split('%').str[0].astype(float) / 100
    # Get depth at each mutation
    s_muts[sp+'_depth'] = s_muts[sp].str.split('(').str[1].str.split(')').str[0].astype(float)
    s_muts[sm+'_depth'] = s_muts[sm].str.split('(').str[1].str.split(')').str[0].astype(float)
    # Get the number of mutant reads for each mutation
    s_muts[sp+'_mutant_reads'] = s_muts[sp+'_af']*s_muts[sp+'_depth']
    s_muts[sm+'_mutant_reads'] = s_muts[sm+'_af']*s_muts[sm+'_depth']
    # Add tumor content as column
    s_muts[sp+'_tc'] = tc.at[sp,'Final tNGS_TC']
    s_muts[sm+'_tc'] = tc.at[sm,'Final tNGS_TC']
    # Calculate expected number of mutant reads per mutation
    s_muts[sp+'_expected_mReads'] = s_muts[sp+'_tc'] * s_muts[sp+'_depth'] / 2
    s_muts[sm+'_expected_mReads'] = s_muts[sm+'_tc'] * s_muts[sm+'_depth'] / 2
    # Add columns for called in each sample
    s_muts[sp+'_called'] = s_muts[sp].str.contains('*', regex = False)
    s_muts[sm+'_called'] = s_muts[sm].str.contains('*', regex = False)
    # set called to true if MAF > 0 and mutant reads > 3
    s_muts.loc[(s_muts[sp+'_af'] > 0) & (s_muts[sp+'_mutant_reads'] >= 3), sp+'_called'] = True
    s_muts.loc[(s_muts[sm+'_af'] > 0) & (s_muts[sm+'_mutant_reads'] >= 3), sm+'_called'] = True
    # Remove rows where expected mutant reads is < 3 in either sample
    s_muts = s_muts[((s_muts[sp+'_expected_mReads'] >= 3)&(s_muts[sm+'_expected_mReads'] >= 3))|((s_muts[sp+'_called'] == True)&(s_muts[sm+'_called'] == True))]
    # Return in no mutations are present
    if len(s_muts) == 0:
        return (np.nan, None)
    # Add shared column. Calculate JSI from this column and return value
    s_muts['shared'] = (s_muts[sp+'_called'] == True)&(s_muts[sm+'_called'] == True)
    s_muts['shared'] = s_muts['shared'].map({True:1,False:0})
    return (s_muts['shared'].sum() / len(s_muts), str(s_muts['shared'].sum())+'/'+str(len(s_muts)))

def get_m2p_jsi(sp, sm, muts):
    s_muts = muts[index_cols+[sp,sm]].copy()
    s_muts = s_muts[(s_muts[sm].str.contains("*", regex = False))]
    # Get the mutant allele frequency of each mutation
    s_muts[sp+'_af'] = s_muts[sp].str.split('%').str[0].astype(float) / 100
    s_muts[sm+'_af'] = s_muts[sm].str.split('%').str[0].astype(float) / 100
    # Get depth at each mutation
    s_muts[sp+'_depth'] = s_muts[sp].str.split('(').str[1].str.split(')').str[0].astype(float)
    s_muts[sm+'_depth'] = s_muts[sm].str.split('(').str[1].str.split(')').str[0].astype(float)
    # Get the number of mutant reads for each mutation
    s_muts[sp+'_mutant_reads'] = s_muts[sp+'_af']*s_muts[sp+'_depth']
    s_muts[sm+'_mutant_reads'] = s_muts[sm+'_af']*s_muts[sm+'_depth']
    # Add tumor content as column
    s_muts[sp+'_tc'] = tc.at[sp,'Final tNGS_TC']
    s_muts[sm+'_tc'] = tc.at[sm,'Final tNGS_TC']
    # Calculate expected number of mutant reads per mutation
    s_muts[sp+'_expected_mReads'] = s_muts[sp+'_tc'] * s_muts[sp+'_depth'] / 2
    s_muts[sm+'_expected_mReads'] = s_muts[sm+'_tc'] * s_muts[sm+'_depth'] / 2
    # Add columns for called in each sample
    s_muts[sp+'_called'] = s_muts[sp].str.contains('*', regex = False)
    s_muts[sm+'_called'] = s_muts[sm].str.contains('*', regex = False)
    # set called to true if MAF > 0 and mutant reads > 3
    s_muts.loc[(s_muts[sp+'_af'] > 0.005) & (s_muts[sp+'_mutant_reads'] >= 3), sp+'_called'] = True
    s_muts.loc[(s_muts[sm+'_af'] > 0.005) & (s_muts[sm+'_mutant_reads'] >= 3), sm+'_called'] = True
    # Remove rows where expected mutant reads is < 3 in either sample
    s_muts = s_muts[((s_muts[sp+'_expected_mReads'] >= 3)&(s_muts[sm+'_expected_mReads'] >= 3))|((s_muts[sp+'_called'] == True)&(s_muts[sm+'_called'] == True))]
    # Return in no mutations are present
    if len(s_muts) == 0:
        return (np.nan, None)
    # Add shared column. Calculate JSI from this column and return value
    s_muts['shared'] = (s_muts[sp+'_called'] == True)&(s_muts[sm+'_called'] == True)
    s_muts['shared'] = s_muts['shared'].map({True:1,False:0})
    return (s_muts['shared'].sum() / len(s_muts), str(s_muts['shared'].sum())+'/'+str(len(s_muts)))

def func(s1, s2):
    if s1.split('_')[2].strip(string.digits) not in s_type_dict.keys():
        print(s1)
    if s2.split('_')[2].strip(string.digits) not in s_type_dict.keys():
        print(s2)
    # Only perform for same patient samples where sample types differ
    if (s1.split('_')[1] != s2.split('_')[1]) or (s_type_dict.get(s1.split('_')[2].strip(string.digits)) == s_type_dict.get(s2.split('_')[2].strip(string.digits))):
        return ((None,None), (None,None));
    
    sp = s1 if s_type_dict.get(s1.split('_')[2].strip(string.digits)) == 'Primary' else (s2 if s_type_dict.get(s2.split('_')[2].strip(string.digits)) == 'Primary' else None)
    sm = s1 if s_type_dict.get(s1.split('_')[2].strip(string.digits)) == 'Metastatic' else (s2 if s_type_dict.get(s2.split('_')[2].strip(string.digits)) == 'Metastatic' else None)
    if pd.isna(sm) or pd.isna(sp):
        return ((None,None), (None,None))
    return get_p2m_jsi(sp, sm, muts), get_m2p_jsi(sp, sm, muts)

m2p_matrix = pd.DataFrame(index = samples, columns = samples)
m2p_matrix_n = pd.DataFrame(index = samples, columns = samples)

p2m_matrix = pd.DataFrame(index = samples, columns = samples)
p2m_matrix_n = pd.DataFrame(index = samples, columns = samples)

s_type_dict = {'RP':'Primary',
               'PB':'Primary',
               'MLN':'Metastatic',
               'MB':'Metastatic',
               'MT':'Metastatic',
               'UCC':'UCC',
               'cfDNA':'Metastatic'}


results = [func(s1,s2) for (s1,s2) in it.combinations_with_replacement(samples, 2)]

n_cores = mp.cpu_count() 
# results = Parallel(n_jobs=n_cores)(delayed(func)(s1,s2) for (s1,s2) in it.combinations_with_replacement(samples, 2))

for (s1,s2),r in zip(it.combinations_with_replacement(samples, 2), results):
    p2m_matrix.at[s1,s2] = r[0][0]
    p2m_matrix.at[s2,s1] = r[0][0]
    p2m_matrix_n.at[s1,s2] = r[0][1]
    p2m_matrix_n.at[s2,s1] = r[0][1]
    
    m2p_matrix.at[s1,s2] = r[1][0]
    m2p_matrix.at[s2,s1] = r[1][0]
    m2p_matrix_n.at[s1,s2] = r[1][1]
    m2p_matrix_n.at[s2,s1] = r[1][1]
    
p2m_matrix.to_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_p2m.tsv', sep = '\t')
p2m_matrix_n.to_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_nMuts_p2m.tsv', sep = '\t')


m2p_matrix.to_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_m2p.tsv', sep = '\t')
m2p_matrix_n.to_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_nMuts_m2p.tsv', sep = '\t')    