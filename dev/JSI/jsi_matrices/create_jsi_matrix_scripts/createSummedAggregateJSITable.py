# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 15:06:55 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np

cohort = 'M1RP'
shared_cols = ['CHROM','POSITION','REF','ALT','GENE','EFFECT']

# =============================================================================
# Import mutations
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
tc = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/tumor_fraction/%s_tumor_fraction_final.tsv' % cohort, sep = '\t')

tc = tc[tc['Final tNGS_TC'] > 0]

# =============================================================================
# Create all sample matrix
# =============================================================================

all_samples = tc['Sample ID'].tolist()
matrix = pd.DataFrame(columns = ['JSI'], index = all_samples)

# =============================================================================
# Select patient and limit mutations to patient
# =============================================================================

for pt in tc['Patient ID'].unique():
    pt_muts = muts[muts['Patient ID'] == pt].copy()
    
    # =============================================================================
    # Create aggregate
    # =============================================================================
    
    agg = pt_muts.copy()
    agg = agg[shared_cols]
    
    for sample in pt_muts['Sample ID'].unique():
        sample_muts = pt_muts[pt_muts['Sample ID'] == sample]
        sample_muts = sample_muts[shared_cols]
        sample_muts['Present'] = 1
        agg = agg.merge(sample_muts, on = shared_cols, how = 'left')
        agg['Present'] = agg['Present'].fillna(0)
        s_jsi = agg['Present'].sum() / len(agg)
        matrix.at[sample, 'JSI'] = s_jsi
        agg = agg.drop('Present', axis = 1)
    
