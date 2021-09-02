# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:09:31 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

    
# =============================================================================
# Import TC
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

# =============================================================================
# Import and melt JSI
# =============================================================================

p2m = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_nMuts_p2m.tsv', sep = '\t', index_col = 'Unnamed: 0')

p2m = p2m.where(np.triu(np.ones(p2m.shape)).astype(bool))
p2m = p2m.reset_index().rename(columns = {'index':'s1'})

p2m = p2m.melt(id_vars = 's1', var_name = 's2', value_name = 'p2m')

p2m = p2m[p2m['s1'] != p2m['s2']]

# =============================================================================
# Import and melt CNV
# =============================================================================

m2p = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES_nMuts_m2p.tsv', sep = '\t', index_col = 'Unnamed: 0')

m2p = m2p.where(np.triu(np.ones(m2p.shape)).astype(bool))
m2p = m2p.reset_index().rename(columns = {'index':'s1'})

m2p = m2p.melt(id_vars = 's1', var_name = 's2', value_name = 'm2p')

m2p = m2p[m2p['s1'] != m2p['s2']]

# =============================================================================
# Merge
# =============================================================================

sim = p2m.merge(m2p, on = ['s1','s2'], how = 'left')

sim['p1'] = sim['s1'].str.split('_').str[1]
sim['p2'] = sim['s2'].str.split('_').str[1]

sim = sim[(~sim['p2m'].isna())&(~sim['m2p'].isna())]

sim = sim[~sim['p1'].isin(['ID8','ID9','ID20'])]
sim = sim[(sim['s1'] != 'M1RP_ID19_cfDNA_2017Jan13')&(sim['s2'] != 'M1RP_ID19_cfDNA_2017Jan13')]

# =============================================================================
# 
# =============================================================================

sim['primary exclusive'] = sim['p2m'].str.split('/').str[1].astype(int) - sim['p2m'].str.split('/').str[0].astype(int) 
sim['metastatic exclusive'] = sim['m2p'].str.split('/').str[1].astype(int) - sim['m2p'].str.split('/').str[0].astype(int) 

# =============================================================================
# 
# =============================================================================

fig,ax1 = plt.subplots(figsize = (2,2))

# =============================================================================
# Plot scatter
# =============================================================================

ax1.scatter(sim['primary exclusive'], sim['metastatic exclusive'], s = 4, lw = 0, c = 'k', zorder = 100)

xlim = ax1.get_xlim()
ylim = ax1.get_ylim()

# ax1.set_xlabel(r'$\frac{\#\ Shared\ mutations}{Primary\ mutations}$')
# ax1.set_ylabel(r'$\frac{\#\ Shared\ mutations}{Metastatic\ mutations}$')

ax1.set_xlabel(r'Number of primary exclusive mutations', fontsize = 6)
ax1.set_ylabel(r'Number of metastatic exclusive mutations', fontsize = 6)

ax1.tick_params(labelsize = 6)

ax1.plot([-2,1000],[-2,1000], ls = 'dashed', zorder = 0, lw = 0.8, c = 'grey')

ax1.set_xlim(xlim)
ax1.set_ylim(ylim)

fig.tight_layout()

ax1.text(100,0,'Wilcoxon p=%.3f' % stats.wilcoxon(sim['primary exclusive'], sim['metastatic exclusive'], alternative='less')[1], fontsize = 6, ha = 'right')

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/exclusive_mutations_priVSmet.pdf')
fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/exclusive_mutations_priVSmet.png')