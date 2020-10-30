# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 14:01:28 2020

@author: amurtha
"""

import pandas as pd
import numpy as np

def sdi(data):
    """ Given a hash { 'species': count } , returns the SDI
    
    >>> sdi({'a': 10, 'b': 20, 'c': 30,})
    1.0114042647073518"""
    
    from math import log as ln
    
    def p(n, N):
        """ Relative abundance """
        if n ==  0:
            return 0
        else:
            return (float(n)/N) * ln(float(n)/N)
            
    N = sum(data.values())
    
    return -sum(p(n, N) for n in data.values() if n != 0)

# =============================================================================
# Import mutation, CN, and TC data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100
tc = tc[tc['Cohort'] == 'M1RP']
tc = tc[tc['Final tNGS_TC'] > 0]

# =============================================================================
# Get patient and sample data. Keep only samples with tumor content
# =============================================================================

all_sdi = {}

for pt in tc['Patient ID'].unique().tolist():
    tc_pt = tc[tc['Patient ID'] == pt].copy()
    muts_pt = muts[muts['Patient ID'] == pt].copy()
    muts_pt = muts_pt[muts_pt['Sample ID'].isin(tc_pt['Sample ID'])]
    
    if len(tc) == 0 or pt == 'ID20': continue;
    
    # =============================================================================
    # Get sets of mutations
    # =============================================================================
    
    unique_muts = muts_pt.drop_duplicates(['CHROM','POSITION','EFFECT']).set_index(['CHROM','POSITION','EFFECT'])
    unique_muts = unique_muts[[]]
    unique_muts['num'] = np.arange(len(unique_muts))
    unique_muts['indentifier'] = 2**unique_muts['num']
    unique_muts = unique_muts.drop('num', axis = 1).reset_index()
    
    muts_pt = muts_pt.merge(unique_muts, on = ['CHROM','POSITION','EFFECT'], how = 'left')
    
    clone_frequencies = muts_pt[['Sample ID','indentifier','Allele_frequency']].groupby('Sample ID').sum().groupby('indentifier').count()
    
    clone_frequencies.columns = ['frequency']
    clone_frequencies['fraction'] = clone_frequencies['frequency'] / clone_frequencies['frequency'].sum()
    
    shannon = sdi(dict(zip(np.arange(len(clone_frequencies)),
                           clone_frequencies['fraction'])))
    all_sdi[pt] = (shannon, clone_frequencies)
    
tmp = pd.DataFrame(all_sdi, index = ['sdi','tmp']).transpose()
tmp = tmp.drop('tmp', axis = 1)
tmp.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/Shannon_index/scripts/targeted_sdi.tsv', sep = '\t')