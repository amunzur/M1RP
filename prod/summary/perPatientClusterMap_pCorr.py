# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 17:13:50 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# =============================================================================
# Import data
# =============================================================================

matrix = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_wesCNA_pearson_corr.tsv', sep = '\t')

melted = matrix.melt(id_vars = ['Unnamed: 0'], var_name = 's2', value_name = 'r')
melted.columns = ['s1','s2','Pearson r']

melted['p1'] = melted['s1'].str.split('_').str[1]
melted['p2'] = melted['s2'].str.split('_').str[1]

melted = melted.groupby(['p1','p2']).median().reset_index()

matrix = melted.pivot(index = 'p1', columns = 'p2', values = 'Pearson r')

grid = sns.clustermap(matrix, figsize=(4, 4), cmap="magma",
                   square = True,
                   col_cluster=True,
                   row_cluster=True,
                   vmin=0,
                   vmax=1,
                   yticklabels=1,
                   xticklabels=1)

grid.ax_heatmap.tick_params(axis = 'both', which = "both", pad=1,
            length = 0, bottom = False, left = False, labelsize=6,
            labelleft = False, labelbottom = False, reset = False)

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/WES_CNA_pCorr_clusterMap_patients.pdf')

fig,ax = plt.subplots()

same = melted[melted['p1'] == melted['p2']].copy()
diff = melted[melted['p1'] != melted['p2']].copy()

diff = diff.drop_duplicates('Pearson r')

same['x'] = same['p1'].apply(lambda x: 1+np.random.uniform(-0.05,0.05))
diff['x'] = diff['p1'].apply(lambda x: 2+np.random.uniform(-0.05,0.05))

ax.boxplot([same['Pearson r'],diff['Pearson r']], showfliers = False, zorder = 100)

melted = pd.concat([same,diff], ignore_index = True)

ax.scatter(melted['x'], melted['Pearson r'], lw = 0, c = 'k', alpha = 0.2, zorder = 10)
