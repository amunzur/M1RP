# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 16:46:40 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# 
# =============================================================================

matrix = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_wesCNA_shared_matrix.tsv', sep = '\t').rename(columns = {'Unnamed: 0':'s1'})

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

sim = matrix.melt(id_vars = 's1', var_name = 's2', value_name = 'similarity')

sim = sim[sim['s1'] != sim['s2']]

sim['p1'] = sim['s1'].str.split('_').str[1]
sim['p2'] = sim['s2'].str.split('_').str[1]

sim = sim.merge(tc[['Sample ID','Final tNGS_TC']], left_on = 's1', right_on = 'Sample ID')
sim = sim.merge(tc[['Sample ID','Final tNGS_TC']], left_on = 's2', right_on = 'Sample ID')

sim = sim[(sim['Final tNGS_TC_x'] > 0.2)&(sim['Final tNGS_TC_y'] > 0.2)]

# =============================================================================
# Get boxplots
# =============================================================================

same = sim[sim['p1'] == sim['p2']]
diff = sim[sim['p1'] != sim['p2']]

same['x'] = same['p1'].apply(lambda x: 1+np.random.uniform(-0.1,0.1))
diff['x'] = diff['p1'].apply(lambda x: 2+np.random.uniform(-0.1,0.1))

fig,ax = plt.subplots()

ax.boxplot([same['similarity'],diff['similarity']], showfliers = False, zorder = 0)
ax.scatter(same['x'], same['similarity'], lw = 0, alpha = 0.3, s = 5, c = 'k', zorder = 100)
ax.scatter(diff['x'], diff['similarity'], lw = 0, alpha = 0.3, s = 5, c = 'k', zorder = 100)
