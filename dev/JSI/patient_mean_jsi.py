# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 17:23:18 2020

@author: amurtha
"""

import pandas as pd
import numpy as np

jsi = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/JSI/jsi_matrices/jsi_matrix_mutationOnly.tsv', sep = '\t', index_col = 'Unnamed: 0')

# =============================================================================
# Get list of patients
# =============================================================================

pts = jsi.index.str.split('_').str[1].unique().tolist()
pts = [pt for pt in pts if pt != 'ID20']

# =============================================================================
# Calculate mean-patient JSI
# =============================================================================

all_mean_jsi = {}

for pt in pts:
    samples = jsi[jsi.index.str.contains(pt+'_')].index.tolist()
    n = len(samples)
    s = 0
    for i, s1 in enumerate(samples[:-1]):
        for s2 in samples[i+1:]:
            s = s + jsi.at[s1, s2]
    s = 2 * s / (n * (n - 1))
    all_mean_jsi[pt] = (s,n)
    
tmp = pd.DataFrame(all_mean_jsi, index = ['JSI','n samples']).transpose()
tmp.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/JSI/patient_mean_jsi.tsv', sep = '\t')