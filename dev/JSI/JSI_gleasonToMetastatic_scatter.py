# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 14:26:01 2020

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
gleason = gleason[gleason['Sample type'] == 'FFPE tissue']

# =============================================================================
# Merge tumor content, keep tumor content > 40
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float)

tc = tc[tc['Final tNGS_TC'] >= 0]

gleason = gleason.merge(tc[['Sample ID', 'Final tNGS_TC']], on = 'Sample ID')

# =============================================================================
# Keep rows as primary samples and columns as metastatic samples in the JSi matrix
# =============================================================================

metastatic = jsi[(jsi.index.str.contains('MLN'))|(jsi.index.str.contains('MB'))].index.tolist()
primary = jsi[(~jsi.index.str.contains('MLN'))&(~jsi.index.str.contains('MB'))].index.tolist()

jsi = jsi[jsi.index.isin(primary)][metastatic]

# =============================================================================
# Merge gleason onto JSI
# =============================================================================

jsi = jsi.merge(gleason[['Sample ID','GGG']], left_index = True, right_on = 'Sample ID')
jsi['GGG'] = jsi['GGG'].astype(int)

# =============================================================================
# Map colors onto gleason scores
# =============================================================================

jsi['Gleason_color'] = jsi['GGG'].apply(lambda x: gleason_colors.get(x))
jsi['Sample color'] = None
jsi.loc[jsi['Sample ID'].str.contains('RP'), 'Sample color'] = 'Red'
jsi.loc[jsi['Sample ID'].str.contains('PB'), 'Sample color'] = 'Blue'


# =============================================================================
# Remove ID8
# =============================================================================

jsi = jsi[~jsi['Sample ID'].str.contains('ID8_')]

# =============================================================================
# set up plot
# =============================================================================

fig,ax = plt.subplots()

# =============================================================================
# Iterate over columns and print non-null jsi
# =============================================================================

for col in jsi.columns:
    if col in ['GGG','Gleason_color','Sample ID','GGG_jitter','Sample color']: continue;
    tmp = jsi[~jsi[col].isnull()].copy()
    tmp['GGG_jitter'] = tmp['GGG'].apply(lambda x: x + (np.random.random()-0.5)/3)
    ax.scatter(tmp['GGG_jitter'],tmp[col], alpha = 0.5, color = tmp['Sample color'], clip_on = False)
    
# =============================================================================
# Create boxplot
# =============================================================================
    
boxplot = []
for ggg in jsi['GGG'].sort_values().unique().tolist():
    ggg_list = []
    for col in jsi.columns:
        if col in ['GGG','Gleason_color','Sample ID','GGG_jitter', 'Sample color']: continue;
        tmp = jsi[(~jsi[col].isnull())&(jsi['GGG'] == ggg)]
        ggg_list = ggg_list + tmp[col].tolist()
    boxplot.append(ggg_list)
    
ax.boxplot(boxplot)
        
kw = stats.kruskal(boxplot[0],boxplot[1],boxplot[2], boxplot[3], boxplot[4])

ax.text(0.75,0.95, 'Kruskal-Wallis test:\np = %.7f' % (round(kw[1], 4)))
    
# # =============================================================================
# # Aethetics
# # =============================================================================

ax.set_ylabel('JSI (primary -> metastatic)')
ax.set_ylim(-0.02, 1)

ax.set_xlabel('Gleason grade group\n(n comparisons)')
ax.set_xticklabels(['1\n(n=%i)' % len(boxplot[0]),'2\n(n=%i)' % len(boxplot[1]),'3\n(n=%i)' % len(boxplot[2]),'4\n(n=%i)' % len(boxplot[3]),'5\n(n=%i)' % len(boxplot[4])])
ax.set_xlim(0.6,5.4)

plt.tight_layout()

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/heterogeniety indexes/JSI_gleasonToMetastatic_scatter.pdf')