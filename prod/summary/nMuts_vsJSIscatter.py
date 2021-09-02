# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:13:27 2021

@author: amurtha
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# =============================================================================
# Import jsi and nmuts
# =============================================================================

jsi = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES.tsv', sep = '\t', index_col = 'Unnamed: 0')

jsi = jsi.where(np.triu(np.ones(jsi.shape)).astype(bool))
jsi = jsi.reset_index().rename(columns = {'index':'s1'})

jsi = jsi.melt(id_vars = 's1', var_name = 's2', value_name = 'JSI')

jsi = jsi[jsi['s1'] != jsi['s2']]
jsi = jsi[~jsi['JSI'].isna()]

# =============================================================================
# Import n muts
# =============================================================================

nmuts = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_nMuts.tsv', sep = '\t', index_col = 'Unnamed: 0')

nmuts = nmuts.where(np.triu(np.ones(nmuts.shape)).astype(bool))
nmuts = nmuts.reset_index().rename(columns = {'index':'s1'})

nmuts = nmuts.melt(id_vars = 's1', var_name = 's2', value_name = 'nmuts')

nmuts = nmuts[nmuts['s1'] != nmuts['s2']]
nmuts = nmuts[~nmuts['nmuts'].isna()]

# =============================================================================
# Merge 2. Keep same patient
# =============================================================================

jsi = jsi.merge(nmuts, on =['s1','s2'])
jsi['p1'] = jsi['s1'].str.split('_').str[1]
jsi['p2'] = jsi['s2'].str.split('_').str[1]

jsi = jsi[jsi['p1'] == jsi['p2']]
jsi = jsi[~jsi['p1'].isin(['ID8','ID9','ID20'])]

jsi = jsi[(jsi['s1'] != 'M1RP_ID19_cfDNA_2017Jan13')&(jsi['s2'] != 'M1RP_ID19_cfDNA_2017Jan13')]
jsi = jsi[(jsi['s1'] != 'M1RP_ID30_UCC')&(jsi['s2'] != 'M1RP_ID30_UCC')]

# =============================================================================
# Plot scatter
# =============================================================================

fig,ax = plt.subplots(figsize = (1.5,2))

ax.scatter(jsi['nmuts'], jsi['JSI'], s = 8, alpha = 0.6, lw = 0, color = 'k')

lin = stats.linregress(jsi['nmuts'], jsi['JSI'])

xlim = ax.get_xlim()
ylim = ax.get_ylim()

b = lin[1]
m = lin[0]

# ax.plot([-10,1000], [b+m*-10, b+m*1000], ls = 'dashed', lw = 0.8, color = 'k', marker = None, zorder)
ax.text(xlim[1], 1.0, '$r^{2}=%0.2f$' % (lin[2]**2), fontsize = 6, ha = 'right')

ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax.set_xlabel('Number of mutations', fontsize = 6)
ax.set_ylabel('JSI', fontsize = 6)
ax.tick_params(labelsize = 6)

plt.tight_layout()
fig.subplots_adjust(left = 0.27, bottom = 0.2)

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/nMuts_vs_JSI.pdf')
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/nMuts_vs_JSI.png')
