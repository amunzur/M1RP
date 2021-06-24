# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:06:39 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import itertools as it
import numpy as np
import multiprocessing as mp
from joblib import Parallel, delayed

effect_dict = {"3'-UTR":'Non-coding', 
               "5'-UTR":'Non-coding', 
               'Downstream':'Non-coding', 
               'Exonic':'Non-coding', 
               'Frameshift':'Truncating', 
               'Intergenic':'Non-coding',
               'Intronic':'Non-coding', 
               'Missense':'Missense', 
               'Non-frameshift':'Non-FS indel', 
               'Splice':'Truncating', 
               'Startloss':'Truncating',
               'Stopgain':'Truncating', 
               'Synonymous':'Non-coding',
               'Upstream':'Non-coding', 
               'Upstream.':'Non-coding'}

index_cols = ['CHROM', 'POSITION', 'GENE', 'REF','ALT', 'EFFECT','NOTES']

def get_jsi(s1, s2, muts):
    s_muts = muts[index_cols+[s1,s2]].copy()
    s_muts = s_muts[(s_muts[s1].str.contains("*", regex = False))|(s_muts[s2].str.contains("*", regex = False))]
    s_muts['s1'] = s1
    s_muts['s2'] = s2
    # Get the mutant allele frequency of each mutation
    s_muts['s1_af'] = s_muts[s1].str.split('%').str[0].astype(float) / 100
    s_muts['s2_af'] = s_muts[s2].str.split('%').str[0].astype(float) / 100
    # Get depth at each mutation
    s_muts['s1_depth'] = s_muts[s1].str.split('(').str[1].str.split(')').str[0].astype(float)
    s_muts['s2_depth'] = s_muts[s2].str.split('(').str[1].str.split(')').str[0].astype(float)
    # Get the number of mutant reads for each mutation
    s_muts['s1_mutant_reads'] = s_muts['s1_af']*s_muts['s1_depth']
    s_muts['s2_mutant_reads'] = s_muts['s2_af']*s_muts['s2_depth']
    # Add tumor content as column
    s_muts['s1_tc'] = tc.at[s1,'Final tNGS_TC']
    s_muts['s2_tc'] = tc.at[s2,'Final tNGS_TC']
    # Calculate expected number of mutant reads per mutation
    s_muts['s1_expected_mReads'] = s_muts['s1_tc'] * s_muts['s1_depth'] / 2
    s_muts['s2_expected_mReads'] = s_muts['s2_tc'] * s_muts['s2_depth'] / 2
    # Add columns for called in each sample
    s_muts['s1_called'] = s_muts[s1].str.contains('*', regex = False)
    s_muts['s2_called'] = s_muts[s2].str.contains('*', regex = False)
    # set called to true if MAF > 0 and mutant reads > 3
    s_muts.loc[(s_muts['s1_af'] > 0) & (s_muts['s1_mutant_reads'] >= 3), 's1_called'] = True
    s_muts.loc[(s_muts['s2_af'] > 0) & (s_muts['s2_mutant_reads'] >= 3), 's2_called'] = True
    # Remove rows where expected mutant reads is < 3 in either sample
    # s_muts = s_muts[(s_muts['s1_expected_mReads'] >= 3)&(s_muts['s2_expected_mReads'] >= 3)]
    # Return in no mutations are present
    if len(s_muts) == 0:
        return np.nan, np.nan
    s_muts = s_muts.drop([s1,s2], axis = 1)
    # Add shared column. Calculate JSI from this column and return value
    s_muts['shared'] = (s_muts['s1_called'] == True)&(s_muts['s2_called'] == True)
    s_muts['shared'] = s_muts['shared'].map({True:1,False:0})
    unshared = s_muts[s_muts['shared'] == 0].copy()
    shared = s_muts[s_muts['shared'] == 1].copy()
    return unshared, shared
    
    
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

samples = muts.columns.tolist()[7:]
results_unshared = []
results_shared = []

for (s1,s2) in it.combinations_with_replacement(samples, 2):
    if s1.split('_')[1] == s2.split('_')[1] and s1 != s2:
        print(s1,s2)
        r = get_jsi(s1, s2, muts)
        if isinstance(r[0], pd.DataFrame) and isinstance(r[1], pd.DataFrame):
            results_unshared.append(r[0])
            results_shared.append(r[1])
            
            
unshared = pd.concat(results_unshared, ignore_index = True)
shared = pd.concat(results_shared, ignore_index = True)

# =============================================================================
# Add patient id
# =============================================================================

unshared['Patient ID'] = unshared['s1'].str.split('_').str[1]
shared['Patient ID'] = shared['s1'].str.split('_').str[1]

# Remove hypermutant from analysis
unshared = unshared[unshared['Patient ID'] != 'ID8']
shared = shared[shared['Patient ID'] != 'ID8']

# Remove ID19 cfDNA lung sample and UCC sample
unshared = unshared[(unshared['s1'] != 'M1RP_ID19_cfDNA_2017Jan13')&(unshared['s2'] != 'M1RP_ID19_cfDNA_2017Jan13')]
unshared = unshared[(unshared['s1'] != 'M1RP_ID30_UCC')&(unshared['s2'] != 'M1RP_ID30_UCC')]

shared = shared[(shared['s1'] != 'M1RP_ID19_cfDNA_2017Jan13')&(shared['s2'] != 'M1RP_ID19_cfDNA_2017Jan13')]
shared = shared[(shared['s1'] != 'M1RP_ID30_UCC')&(shared['s2'] != 'M1RP_ID30_UCC')]

# =============================================================================
# Remove duplicates. If unshared, it is not shared
# =============================================================================

unshared = unshared.drop_duplicates(['CHROM', 'POSITION', 'GENE', 'REF', 'ALT', 'Patient ID'])
shared = shared.drop_duplicates(['CHROM', 'POSITION', 'GENE', 'REF', 'ALT', 'Patient ID'])

unshared = unshared.set_index(['CHROM', 'POSITION', 'GENE', 'REF', 'ALT', 'Patient ID'])
shared = shared.set_index(['CHROM', 'POSITION', 'GENE', 'REF', 'ALT', 'Patient ID'])

unshared = unshared[(unshared['s1_expected_mReads'] >= 3)&(unshared['s2_expected_mReads'] >= 3)|((unshared['s1_called'] == True)&(unshared['s2_called'] == True))]
shared = shared[(shared['s1_expected_mReads'] >= 3)&(shared['s2_expected_mReads'] >= 3)|((shared['s1_called'] == True)&(shared['s2_called'] == True))]

shared = shared[~shared.index.isin(unshared.index)]

# =============================================================================
# 
# =============================================================================

unshared['EFFECT'] = unshared['EFFECT'].str.split(' ').str[0]
shared['EFFECT'] = shared['EFFECT'].str.split(' ').str[0]

unshared['EFFECT'] = unshared['EFFECT'].map(effect_dict)
shared['EFFECT'] = shared['EFFECT'].map(effect_dict)


mut_type_unshared = unshared.groupby(['EFFECT']).count()
mut_type_shared = shared.groupby(['EFFECT']).count()

mut_type = mut_type_shared[['s1']].merge(mut_type_unshared[['s2']], left_index = True, right_index = True)
mut_type.columns = ['shared_count','unshared_count']

# mut_type['shared_count'] = mut_type['shared_count'] / len(shared)
# mut_type['unshared_count'] = mut_type['unshared_count'] / len(unshared)

# =============================================================================
# Plot
# =============================================================================
    
fig,ax = plt.subplots()

ax.bar(np.arange(0.0,len(mut_type),1),mut_type['shared_count'], width = 0.4, label = 'Truncal', color = 'k')
ax.bar(np.arange(0.4,len(mut_type),1),mut_type['unshared_count'], width = 0.4, label = 'Non-truncal', color = 'grey')

ax.set_xticks(np.arange(0.2,len(mut_type),1))
ax.set_xticklabels(mut_type.index)

ax.legend()



'''
def func(s1, s2):
    if s1.split('_')[1] == s2.split('_')[1] and s1 != s2:
        return get_jsi(s1, s2, muts)
    return

results = [func(s1,s2) for (s1,s2) in it.combinations_with_replacement(samples, 2)]

n_cores = mp.cpu_count() 
results = Parallel(n_jobs=n_cores)(delayed(func)(s1,s2) for (s1,s2) in it.combinations_with_replacement(samples, 2))
'''