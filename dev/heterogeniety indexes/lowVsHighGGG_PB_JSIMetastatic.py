# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 23:41:14 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import rgb2hex
import scipy.stats as stats

cmap=plt.get_cmap('OrRd', 6)
gleason_colors = dict(zip(np.arange(1,6,1), [rgb2hex(cmap(int(i+1))) for i in np.arange(5)]))


# =============================================================================
# Import data. Gleason + JSI
# =============================================================================

gleason = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')
jsi = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/jsi_matrices/jsi_matrix_mutationOnly.tsv', sep = '\t', index_col = 'Unnamed: 0')

gleason.columns = gleason.iloc[0]
gleason = gleason.drop(gleason.index[0])

# =============================================================================
# Keep M1RP tissue samples only
# =============================================================================

gleason = gleason[gleason['Cohort'] == 'M1RP']
gleason = gleason[~gleason['Sample ID'].isnull()]
gleason = gleason[gleason['Sample ID'].str.contains('PB')]

gleason = gleason[gleason['Patient ID'] != 'ID8']

# =============================================================================
# Merge tumor content onto gleason dataframe. Keep TC > 0
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float)

tc = tc[tc['Final tNGS_TC'] >= 0]

gleason = gleason.merge(tc[['Sample ID', 'Final tNGS_TC']], on = 'Sample ID')

# =============================================================================
# Get patients with more than 1 gleason score witin PB samples
# =============================================================================

pt_counts = gleason.groupby('Patient ID').nunique()
pt_counts = pt_counts[pt_counts['GGG'] > 1]
gleason = gleason[gleason['Patient ID'].isin(pt_counts.index.tolist())]

# =============================================================================
# Get lowest and highest gleason scores
# =============================================================================

gleason_low = gleason.sort_values(['Patient ID', 'GGG','Final tNGS_TC'], ascending = True)
gleason_high = gleason.sort_values(['Patient ID', 'GGG','Final tNGS_TC'], ascending = False)

gleason_low = gleason_low.drop_duplicates(['Patient ID'])
gleason_high = gleason_high.drop_duplicates(['Patient ID'])

# =============================================================================
# Keep only metastatic samples in JSI columns with TC > 0
# =============================================================================

met_samples = tc.copy()
met_samples = met_samples[(met_samples['Sample ID'].str.contains('MLN'))|(met_samples['Sample ID'].str.contains('MB'))|(met_samples['Sample ID'].str.contains('cfDNA'))]
met_samples = met_samples[met_samples['Cohort'] == 'M1RP']
met_samples = met_samples[met_samples['Final tNGS_TC'] > 0]
met_samples = met_samples[met_samples['Patient ID'].isin(pt_counts.index.tolist())]

# =============================================================================
# Create low and high JSI matrices
# =============================================================================

jsi['Patient ID'] = jsi.index.str.split('_').str[1]
jsi = jsi[jsi['Patient ID'].isin(pt_counts.index)]
jsi = jsi[met_samples['Sample ID'].tolist()]

jsi_low = jsi[jsi.index.isin(gleason_low['Sample ID'])].copy()
jsi_high = jsi[jsi.index.isin(gleason_high['Sample ID'])].copy()

# =============================================================================
# Remove all NA rows
# =============================================================================

jsi_low = jsi_low.dropna(how = 'all')
jsi_high = jsi_high.dropna(how = 'all')

# =============================================================================
# Take median JSI for each samples
# =============================================================================

jsi_low['Median JSI'] = jsi_low.median(axis = 1, skipna = True)
jsi_high['Median JSI'] = jsi_high.median(axis = 1, skipna = True)

# =============================================================================
# Create plot
# =============================================================================

fig,ax = plt.subplots()

# =============================================================================
# Plot scatter for low and high
# =============================================================================

jitter_low = (np.random.rand(len(jsi_low))/10+0.9)
ax.scatter(jitter_low, jsi_low['Median JSI'], color = 'k')

jitter_high = (np.random.rand(len(jsi_high))/10+1.9)
ax.scatter(jitter_high, jsi_high['Median JSI'], color = 'k')

# =============================================================================
# Plot boxplots
# =============================================================================

ax.boxplot([jsi_low['Median JSI'], jsi_high['Median JSI']])

# =============================================================================
# Aethetics
# =============================================================================

ax.set_xticklabels(['Lowest GGG\nprostate biopsy','Highest GGG\nprostate biopsy'])
ax.set_ylabel('Median metastatic JSI')

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/lowVsHighGGG_PB_JSIMetastatic.pdf')