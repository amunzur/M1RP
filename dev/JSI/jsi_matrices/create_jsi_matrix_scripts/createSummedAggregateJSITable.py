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
        if sample == 'M1RP_ID2_RP3':
            sample = sample
        sample_muts = pt_muts[pt_muts['Sample ID'] == sample]
        sample_muts = sample_muts[shared_cols]
        sample_muts['Present'] = 1
        agg = agg.merge(sample_muts, on = shared_cols, how = 'left')
        agg['Present'] = agg['Present'].fillna(0)
        s_jsi = agg['Present'].sum() / len(agg)
        matrix.at[sample, 'JSI'] = s_jsi
        agg = agg.drop('Present', axis = 1)
    
matrix = matrix.reset_index().rename(columns = {'index':'Sample ID'})
matrix['Patient ID'] = matrix['Sample ID'].str.split('_').str[1]

# Create X coordinates to be plotted
matrix['x'] = matrix['Patient ID'].str.split('ID').str[1].astype(int)
matrix['x'] = matrix['x'].apply(lambda x: x + (random.random()-0.5)/5)

matrix = matrix[['Patient ID','Sample ID','JSI','x']]
matrix = matrix[matrix['Sample ID'].isin(tc['Sample ID'])]

# =============================================================================
# Create boxplot arrays
# =============================================================================
    
bp = []
for pt in matrix['Patient ID'].unique().tolist():
    bp.append(matrix[matrix['Patient ID'] == pt]['JSI'].tolist())

# =============================================================================
# Create swarmplot of JSI values per patient
# =============================================================================

fig,ax = plt.subplots()

# ax.scatter(matrix['x'], matrix['JSI'], c = 'k', alpha = 0.6, s = 8)
ax.boxplot(bp)
ax.set_xticks(np.arange(1,matrix['x'].max()+0.25,1))