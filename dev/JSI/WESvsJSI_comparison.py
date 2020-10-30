# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 15:59:49 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import itertools
import scipy.stats as stats

# =============================================================================
# Import JSI values
# =============================================================================

targeted = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/JSI/jsi_matrices/jsi_matrix_mutThenCN.tsv', sep = '\t', index_col = 'Unnamed: 0')
wes = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/JSI/jsi_matrices/jsi_matrix_wxs.tsv', sep = '\t', index_col = 'Unnamed: 0')

# =============================================================================
# Keep only samples in WES
# =============================================================================

samples = wes.index.tolist()
targeted = targeted[targeted.index.isin(samples)]

samples = targeted.index.tolist()
targeted = targeted[samples]
wes = wes[wes.index.isin(samples)][samples]

# =============================================================================
# Fill out dataframe of combinations
# =============================================================================

comb = pd.DataFrame(index = pd.MultiIndex.from_tuples(list(itertools.combinations(samples ,2))))

for index, row in comb.iterrows():
    s1 = index[0]
    s2 = index[1]
    comb.at[index, 'Targeted'] = targeted.at[s1,s2]
    comb.at[index, 'WES'] = wes.at[s1,s2]
    
comb = comb.dropna(how = 'all')
    
# =============================================================================
# Plot data
# =============================================================================

fig,ax = plt.subplots()

ax.scatter(comb['Targeted'], comb['WES'], s = 7, color = 'k', alpha = 0.5)

ax.set_ylabel('WES JSI')
ax.set_xlabel('Targeted_JSI')

# =============================================================================
# Get the linear stats of data and plot line and r value
# =============================================================================

lin = stats.linregress(comb['Targeted'], comb['WES'])
ax.text(x = 0, y = 0.52, s = 'r: %.2f\np: %.3f' % (lin[2], lin[3]))

# Plot trendline
ax.plot([-1,1], [-1 * lin[0] + lin[1], 1*lin[0]+lin[1]])

ax.set_ylim(-.025, 0.60)
ax.set_xlim(-.025, 0.95)

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/JSI/WESvsJSI_comparison.pdf')