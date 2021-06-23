# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 17:43:02 2021

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

sim = sim[(sim['Final tNGS_TC_x'] > 0.2)&(sim['Final tNGS_TC_y'] > 0.2)]

sim.loc[sim['JSI'] < 0, 'JSI'] = 0

# =============================================================================
# Get boxplots for same patient vs other patient
# =============================================================================

same = sim[sim['p1'] == sim['p2']].copy()
diff = sim[sim['p1'] != sim['p2']].copy()

secondary_tumors = ['M1RP_ID19_cfDNA_2017Jan13','M1RP_ID30_UCC']

same = same[(same['p1']!='ID8')&(same['p2']!='ID8')]
same = same[(~same['s1'].isin(secondary_tumors))&(~same['s2'].isin(secondary_tumors))]

same['x'] = same['p1'].apply(lambda x: 1+np.random.uniform(-0.1,0.1))
diff['x'] = diff['p1'].apply(lambda x: 2+np.random.uniform(-0.1,0.1))

fig,ax = plt.subplots(figsize=(1.5,2))

ax.boxplot([same['JSI'],diff['JSI']], showfliers = False, zorder = 0)
ax.scatter(same['x'], same['JSI'], lw = 0, alpha = 0.3, s = 5, c = 'k', zorder = 100)
ax.scatter(diff['x'], diff['JSI'], lw = 0, alpha = 0.1, s = 5, c = 'k', zorder = 100)

ax.set_xticklabels(['Same\npatient','Different\npatient'])
ax.set_ylabel('WES JSI', fontsize = 6)

ax.tick_params(labelsize = 6)

fig.tight_layout()

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/wes_snv_JSI_boxplots_samePt_diffPt.pdf')


# =============================================================================
# Get JSI between metastatic and primary
# =============================================================================

t_dict = {'RP':'Primary', 'MLN':'Metastatic', 'cfDNA':'cfDNA', 'PB':'Primary', 'MB':'Metastatic','UCC':'UCC'}

same['t1'] = same['s1'].str.split('_').str[2].str.strip(string.digits).map(t_dict)
same['t2'] = same['s2'].str.split('_').str[2].str.strip(string.digits).map(t_dict)

same = same[(same['p1']!='ID8')&(same['p2']!='ID8')]

# =============================================================================
# Plot primary primary, met_primary, met_met, cfDNA_primary, cfDNA_met
# =============================================================================

pp = same[(same['t1'] == 'Primary')&(same['t2'] == 'Primary')].copy()
mp = same[((same['t1'] == 'Primary')&(same['t2'] == 'Metastatic'))|((same['t1'] == 'Metastatic')&(same['t2'] == 'Primary'))].copy()
mm = same[(same['t1'] == 'Metastatic')&(same['t2'] == 'Metastatic')].copy()
cp = same[((same['t1'] == 'Primary')&(same['t2'] == 'cfDNA'))|((same['t1'] == 'cfDNA')&(same['t2'] == 'Primary'))].copy()
mc = same[((same['t1'] == 'cfDNA')&(same['t2'] == 'Metastatic'))|((same['t1'] == 'Metastatic')&(same['t2'] == 'cfDNA'))].copy()

pp['x'] = pp['p1'].apply(lambda x: 1+np.random.uniform(-0.1,0.1))
mp['x'] = mp['p1'].apply(lambda x: 2+np.random.uniform(-0.1,0.1))
mm['x'] = mm['p1'].apply(lambda x: 3+np.random.uniform(-0.1,0.1))
cp['x'] = cp['p1'].apply(lambda x: 4+np.random.uniform(-0.1,0.1))
mc['x'] = mc['p1'].apply(lambda x: 5+np.random.uniform(-0.1,0.1))

# =============================================================================
# Plot
# =============================================================================

fig,ax = plt.subplots(figsize = (3,2))

ax.boxplot([pp['JSI'],mp['JSI'],mm['JSI'],cp['JSI'],mc['JSI']], showfliers = False, zorder = 0)
ax.scatter(pp['x'], pp['JSI'], lw = 0, alpha = 0.5, s = 8, c = 'k', zorder = 100)
ax.scatter(mp['x'], mp['JSI'], lw = 0, alpha = 0.5, s = 8, c = 'k', zorder = 100)
ax.scatter(mm['x'], mm['JSI'], lw = 0, alpha = 0.5, s = 8, c = 'k', zorder = 100)
ax.scatter(cp['x'], cp['JSI'], lw = 0, alpha = 0.5, s = 8, c = 'k', zorder = 100)
ax.scatter(mc['x'], mc['JSI'], lw = 0, alpha = 0.5, s = 8, c = 'k', zorder = 100)

labels = ['Primary\nPrimary\nn=%i' % len(pp),'Primary\nmetastatic\nn=%i' % len(mp),'Metastatic\nmetastatic\nn=%i' % len(mm),'Primary\ncfDNA\nn=%i' % len(cp),'Metastatic\ncfDNA\nn=%i' % len(mc)]

ax.set_xticklabels(labels)
ax.set_ylabel('WES JSI', fontsize = 6)

ax.tick_params(labelsize = 6)

plt.tight_layout()

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/wes_snv_JSI_boxplots_bySite.pdf')
