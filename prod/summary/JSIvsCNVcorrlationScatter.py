# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:21:56 2021

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

jsi = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES.tsv', sep = '\t', index_col = 'Unnamed: 0')

jsi = jsi.where(np.triu(np.ones(jsi.shape)).astype(bool))
jsi = jsi.reset_index().rename(columns = {'index':'s1'})

jsi = jsi.melt(id_vars = 's1', var_name = 's2', value_name = 'JSI')

jsi = jsi[jsi['s1'] != jsi['s2']]

# =============================================================================
# Import and melt CNV
# =============================================================================

cnv = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_wesCNA_pearson_corr.tsv', sep = '\t', index_col = 'Unnamed: 0')

cnv = cnv.where(np.triu(np.ones(cnv.shape)).astype(bool))
cnv = cnv.reset_index().rename(columns = {'index':'s1'})

cnv = cnv.melt(id_vars = 's1', var_name = 's2', value_name = 'p_corr')

cnv = cnv[cnv['s1'] != cnv['s2']]

# =============================================================================
# Get boxplots for same patient vs other patient
# =============================================================================

sim = cnv.merge(jsi, on = ['s1','s2'], how = 'left')

sim['p1'] = sim['s1'].str.split('_').str[1]
sim['p2'] = sim['s2'].str.split('_').str[1]

sim = sim.merge(tc[['Sample ID','Final tNGS_TC']], left_on = 's1', right_on = 'Sample ID')
sim = sim.merge(tc[['Sample ID','Final tNGS_TC']], left_on = 's2', right_on = 'Sample ID')

# =============================================================================
# Keep high tc and same patient combinations
# =============================================================================

sim = sim[(sim['Final tNGS_TC_x'] > 0.2)&(sim['Final tNGS_TC_y'] > 0.2)]
sim = sim[sim['p1'] == sim['p2']]
sim = sim[(~sim['p_corr'].isna())&(~sim['JSI'].isna())]

secondary_tumors = ['M1RP_ID19_cfDNA_2017Jan13','M1RP_ID30_UCC']

sim = sim[sim['p1'] != 'ID8']
sim = sim[(~sim['s1'].isin(secondary_tumors))&(~sim's2'].isin(secondary_tumors))]

# =============================================================================
# Plot scatter
# =============================================================================

fig,ax = plt.subplots()

ax.scatter(sim['JSI'], sim['p_corr'], lw = 0, c = 'k', s = 12)
lin = stats.linregress(sim['JSI'],sim['p_corr'])

xlim = ax.get_xlim()
ylim = ax.get_ylim()

b = lin[1]
m = lin[0]

ax.plot([-10,10], [b+m*-10, b+m*10], ls = 'dashed', lw = 0.8, color = 'k', marker = None)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax.set_xlabel('JSI')
ax.set_ylabel('CNV correlation')