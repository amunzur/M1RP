# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 16:57:59 2021

@author: amurtha
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import string

# =============================================================================
# Import data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_all-segments.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

# =============================================================================
# Select patient and import WES CNA
# =============================================================================

cn = cn.merge(tc[['Sample ID','Final tNGS_TC']], on = 'Sample ID', how = 'left')
cn.loc[cn['Sample ID'] == 'M1RP_ID33_cfDNA_2021Mar30', 'Final tNGS_TC'] = 0.41

cn = cn[cn['Final tNGS_TC'] >= 0.1]

# =============================================================================
# Calculate the estimated copy number
# =============================================================================

cn = cn[cn['CHR'] != 'chrM']
cn = cn[cn['CHR'] != 'chrY']

cn['copies'] = -1
for allosome in [True,False]:
    if allosome:
        cn.loc[cn['CHR'].isin(['chrX',]), 'copies'] = 1+(2**cn['Log_ratio']-1)/cn['Final tNGS_TC']
    else:
        cn.loc[~cn['CHR'].isin(['chrX']), 'copies'] = (2**(cn['Log_ratio']+1)+2*cn['Final tNGS_TC'] - 2) / cn['Final tNGS_TC']

cn.loc[cn['Log_ratio'] == 0, 'copies'] = 2
cn['copies'] = cn['copies'].round()

# =============================================================================
# Calculate the total panel size
# =============================================================================

cn['Patient ID'] = cn['Sample ID'].str.split('_').str[1]
cn['seg_size'] = cn['END'] - cn['START']

pt_size = cn.copy()
pt_size =  pt_size.drop_duplicates(['Sample ID','CHR','START','END'])
pt_size = cn.groupby(['Sample ID']).sum().reset_index()

pt_size = dict(zip(pt_size['Sample ID'],pt_size['seg_size']))

cn['pt_size'] = cn['Sample ID'].map(pt_size)

# =============================================================================
# Calculate % altered
# =============================================================================

cn['altered'] = 1
cn.loc[cn['copies'] == 2, 'altered'] = 0

cn['rel_size'] = cn['seg_size'] / cn['pt_size']
cn['altered_size'] = cn['rel_size']*cn['altered']

cn = cn.groupby('Sample ID').sum().reset_index()

cn_all = cn.copy()

# =============================================================================
# Keep highest TC met and primary and cfDNA samples
# =============================================================================

s_type = {'RP':'Primary', 'MLN':'Metastatic', 'cfDNA':'cfDNA', 'PB':'Primary', 'MB':'Metastatic','UCC':'UCC'}

tc = tc[tc['Sample ID'].isin(cn['Sample ID'].tolist())]

tc['Sample category'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits)
tc['Sample category'] = tc['Sample category'].map(s_type)

tc = tc[tc['Sample ID'] != 'M1RP_ID30_RP3']
tc = tc[tc['Sample category'] != 'UCC']

tc = tc.sort_values('Final tNGS_TC', ascending = False)

tc_m = tc[tc['Sample category'] == 'Metastatic']
tc_p = tc[tc['Sample category'] == 'Primary']
tc_c = tc[tc['Sample category'] == 'cfDNA']

tc_m = tc_m.drop_duplicates('Patient ID')
tc_p = tc_p.drop_duplicates('Patient ID')
tc_c = tc_c.drop_duplicates('Patient ID')

samples = tc_m['Sample ID'].tolist()+tc_p['Sample ID'].tolist()+tc_c['Sample ID'].tolist()+['M1RP_ID33_cfDNA_2021Mar30']

cn = cn[cn['Sample ID'].isin(samples)]

# =============================================================================
# Set up x values
# =============================================================================

cn_all = cn_all.merge(tc[['Sample ID','Sample category']], on = 'Sample ID', how = 'left')
cn = cn.merge(tc[['Sample ID','Sample category']], on = 'Sample ID', how = 'left')
cn.loc[cn['Sample ID'].str.contains('cfDNA'), 'Sample category'] = 'cfDNA'

x_dict = dict(zip(['Primary','Metastatic','cfDNA'],[1,2,3]))

jitter = 0.1
cn['x'] = cn['Sample category'].apply(lambda cat: x_dict.get(cat)+np.random.uniform(low=jitter, high=-1*jitter))

# =============================================================================
# Get boxplots
# =============================================================================

cn_m = cn[cn['Sample category'] == 'Metastatic'].copy()
cn_p = cn[cn['Sample category'] == 'Primary'].copy()
cn_c = cn[cn['Sample category'] == 'cfDNA'].copy()

boxplots = []

for df in [cn_p, cn_m, cn_c]:
    boxplots.append(df['altered_size'].tolist())
    
# =============================================================================
# Plot boxplots
# =============================================================================

fig,[ax1,ax2] = plt.subplots(ncols = 2, figsize = (3.5,2))

ax1.boxplot(boxplots, showfliers = False, zorder = 0)
ax1.scatter(cn['x'],cn['altered_size'], lw = 0, s = 6, zorder = 100,clip_on = False)

# =============================================================================
# Set up the lines connecting same patient
# =============================================================================

cn_m['Patient ID'] = cn_m['Sample ID'].str.split('_').str[1]
cn_p['Patient ID'] = cn_p['Sample ID'].str.split('_').str[1]
cn_c['Patient ID'] = cn_c['Sample ID'].str.split('_').str[1]

cn_delta = pd.DataFrame({'Patient ID':cn_p['Patient ID'].unique()})
cn_delta = cn_delta.merge(cn_p[['Patient ID','altered_size','x']], on = 'Patient ID', how = 'left')
cn_delta = cn_delta.merge(cn_m[['Patient ID','altered_size','x']], on = 'Patient ID', how = 'left')
cn_delta = cn_delta.merge(cn_c[['Patient ID','altered_size','x']], on = 'Patient ID', how = 'left')

cn_delta.columns = ['Patient ID','Primary','p_x','Metastatic','m_x','cfDNA','c_x',]

p_m = cn_delta[~cn_delta['Metastatic'].isna()]
m_c = cn_delta[(~cn_delta['Metastatic'].isna())&(~cn_delta['cfDNA'].isna())]

ax2.scatter(p_m['Primary'],p_m['Metastatic'], c = 'k', zorder = 100, lw = 0, s = 8, clip_on = False)
ax2.plot([0,1],[0,1], lw = 0.8, ls = 'dashed', zorder = 0, color = 'k', alpha = 0.4)

# =============================================================================
# Aethetics
# =============================================================================

# AX1
ax1.set_ylim(0,0.6)
ax1.set_ylabel('Percent genome altered', fontsize = 6)

ns = [len(x) for x in boxplots]
ax1.set_xticklabels(['Primary\nn=%i' % ns[0],'Metastatic\nn=%i' % ns[1],'cfDNA\nn=%i' % ns[2]])

ax1.tick_params(labelsize = 6)

# Ax2
ax2.set_ylim(0,0.6)
ax2.set_ylabel('Percent genome altered (Metastatic)', fontsize = 6)

ax2.set_xlim(0,0.6)
ax2.set_xlabel('Percent genome altered (Priamry)', fontsize = 6)

ax2.tick_params(labelsize = 6)
fig.tight_layout()

import scipy.stats as stats

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/WES_percent_altered.pdf')

cn_all = cn_all[['Sample ID','altered_size','Sample category']]
cn_all.columns = ['Sample ID','Percent altered','Sample category']

cn_all.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/M1RP Copy Number Analysis/Whole exome (data files, threshold 2.0)/M1RP_percentAltered.tsv')

st = stats.ttest_rel(p_m['Primary'],p_m['Metastatic'], alternative = 'less')

