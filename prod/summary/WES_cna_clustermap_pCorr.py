# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 14:14:47 2021

@author: amurtha
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# =============================================================================
# 
# =============================================================================

matrix = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_wesCNA_pearson_corr.tsv', sep = '\t', index_col = 'Unnamed: 0')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_all-segments.tsv', sep = '\t')

# =============================================================================
# Keep samples from patients with more than 1 sample
# =============================================================================

samples = matrix.index.tolist()

df = pd.DataFrame({'Sample ID':samples})

df['Patient ID'] = df['Sample ID'].str.split('_').str[1]

g = df.groupby(['Patient ID']).count()
g = g[g['Sample ID'] > 1]
df = df[df['Patient ID'].isin(g.index.tolist())]

matrix = matrix[matrix.index.isin(df['Sample ID'].tolist())]
matrix = matrix[df['Sample ID'].tolist()]

# =============================================================================
# Remove samples with 0 copy number alterations
# =============================================================================

# mask = matrix.isnull().copy()

# matrix = matrix.fillna(0)

# =============================================================================
# Create heatmap
# =============================================================================

grid = sns.clustermap(matrix, figsize=(14, 14), cmap="magma",
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

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/WES_CNA_pCorr_clusterMap.pdf')