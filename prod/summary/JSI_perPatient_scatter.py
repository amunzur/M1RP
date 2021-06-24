# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 17:11:04 2021

@author: amurtha
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string

# =============================================================================
# 
# =============================================================================

matrix = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/prod/heterogeniety indexes/jsi_matrix_WES.tsv', sep = '\t', index_col = 'Unnamed: 0')

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

matrix = matrix.where(np.triu(np.ones(matrix.shape)).astype(bool))
matrix = matrix.reset_index().rename(columns = {'index':'s1'})

sim = matrix.melt(id_vars = 's1', var_name = 's2', value_name = 'JSI')

sim = sim[~sim['JSI'].isna()]
sim = sim[sim['s1'] != sim['s2']]

sim['p1'] = sim['s1'].str.split('_').str[1]
sim['p2'] = sim['s2'].str.split('_').str[1]

sim = sim.merge(tc[['Sample ID','Final tNGS_TC']], left_on = 's1', right_on = 'Sample ID')
sim = sim.merge(tc[['Sample ID','Final tNGS_TC']], left_on = 's2', right_on = 'Sample ID')

sim = sim[(sim['Final tNGS_TC_x'] > 0.1)&(sim['Final tNGS_TC_y'] > 0.1)]

# =============================================================================
# Get boxplots for same patient vs other patient
# =============================================================================

sim = sim[sim['p1'] == sim['p2']].copy()

secondary_tumors = ['M1RP_ID19_cfDNA_2017Jan13','M1RP_ID30_UCC']

sim['color'] = 'k'
sim.loc[sim['p1'] == 'ID8', 'color'] = 'orange'
sim.loc[(sim['s1'].isin(secondary_tumors))|(sim['s2'].isin(secondary_tumors)), 'color'] = 'blue'

# =============================================================================
# Get order
# =============================================================================

pt_median = sim.groupby('p1').median().sort_values('JSI', ascending = False)
x_dict = dict(zip(pt_median.index.tolist(), np.arange(len(pt_median))))

sim['x'] = sim['p1'].map(x_dict)

# =============================================================================
# plot
# =============================================================================

fig,ax = plt.subplots()

ax.scatter(sim['x'],sim['JSI'], alpha = 0.6, lw = 0, s = 11, color = sim['color'])

ax.set_xticks(np.arange(len(pt_median)))
ax.set_xticklabels(list(x_dict.keys()), rotation = 90)

ax.set_ylabel('JSI')