# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 12:15:42 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import random
import numpy as np

'''
Figure will be a scatter plot of mean CNA deviation from 0 vs TC where TC is a 
rolling average 
'''

# =============================================================================
# Helpers
# =============================================================================

# Return theoretical log-ratio for a given TC and copy number
def getLRFromCNandTC(tc, cn):
    tc = cn * tc - 2 * tc + 2
    tc = np.log2(tc) - 1
    return np.abs(tc)
    

# =============================================================================
# Import data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_FFPE_cna.tsv', sep = '\t')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc['Group'] = tc['Group'].astype(int)

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

# =============================================================================
# Keep only tissue samples
# =============================================================================

samples = samples[samples['Cohort'] == 'M1RP']
samples = samples[samples['Sample Category'].isin(['PB','RP','MLN','MB','MT'])]

tc = tc[tc['Sample ID'].isin(samples['Sample ID'].tolist())]

# =============================================================================
# Keep only samples with > 0 TC
# =============================================================================

tc = tc[tc['Final tNGS_TC'] > 0]

cn = cn.merge(tc[['Sample ID', 'Final tNGS_TC']], on = 'Sample ID')

# =============================================================================
# Get ABS of log-ratio
# =============================================================================

bad_genes = ['CDKN1B','AKT1','ERCC2','AKT3','HSD3B1','TMPRSS2']

# cn = cn[cn['Copy_num'] != 0]
cn = cn[~cn['GENE'].isin(bad_genes)]
cn['Log_ratio'] = cn['Log_ratio'].abs()

# =============================================================================
# Create plot
# =============================================================================

fig,[ax1,ax2] = plt.subplots(ncols = 2)

# =============================================================================
# Create rolling average. Plot scatter
# =============================================================================

xs = []
ys = []

for x in np.arange(0.05, 1.0, 0.05):
    tmp = cn[(cn['Final tNGS_TC'] < x+0.05)&(cn['Final tNGS_TC'] >= x-0.05)].copy()
    
    ax1.scatter(x,tmp['Log_ratio'].mean(), c = 'k')
    # ax1.errorbar(x,tmp['Log_ratio'].mean(), yerr = tmp['Log_ratio'].std()*2)
    xs.append(x)
    ys.append(tmp['Log_ratio'].mean())
    
lin = stats.linregress(xs,ys)

# ax1.set_ylim(0.02, 0.45)

ax1.set_xlabel('Tumor content (rolling average +/-0.05)')
ax1.set_ylabel('Mean absolute log-ratio deviation')

ax1.text(0.01,0.3,'r = %.3f' % lin[2])

# =============================================================================
# Plot boxplots and swarms
# =============================================================================

med_tc = cn.drop_duplicates('Sample ID')['Final tNGS_TC'].median()

above = cn[cn['Final tNGS_TC'] >= med_tc].copy()
below = cn[cn['Final tNGS_TC'] <  med_tc].copy()

above['x'] = 2
above['x'] = above['x'].apply(lambda x: x + random.uniform(-0.4, 0.4))

below['x'] = 1
below['x'] = below['x'].apply(lambda x: x + random.uniform(-0.4, 0.4))

# ax2.boxplot([below['Log_ratio'].tolist(),above['Log_ratio'].tolist()], showfliers = False)
ax2.scatter(above['x'],above['Log_ratio'], s = 2, lw = 0, c = 'k', alpha = 0.3)
ax2.scatter(below['x'],below['Log_ratio'], s = 2, lw = 0, c = 'k', alpha = 0.3)

ax2.text(0.6, 1.3, 'Mann-Whitney U (1 sided)\np = %.3f\nTC median = %.3f' % (stats.mannwhitneyu(above['Log_ratio'], below['Log_ratio'], alternative='greater')[1], med_tc))

ax2.set_xticks([1,2])
ax2.set_xticklabels(['TC < median', 'TC >= median'])
ax2.set_ylim(0,1.5)
ax2.set_ylabel('Mean absolute log-ratio deviation')

fig.savefig('G:/Andy Murtha/Ghent/M1RP/prod/summary/CNlr_vs_TC.pdf')