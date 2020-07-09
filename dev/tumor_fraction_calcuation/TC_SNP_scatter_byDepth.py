# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 14:20:01 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats

# =============================================================================
# import data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
depth = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=878831703')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['snp_TC'] = (tc['snp_TC'].str.split('%').str[0].astype(float) / 100).fillna(0)
tc['mut_TC'] = (tc['mut_TC'].str.split('%').str[0].astype(float) / 100).fillna(0)

# =============================================================================
# Keep only targeted non gDNA data
# =============================================================================

depth = depth[depth['Targeted or WES'] == 'Targeted (73 gene)']
depth = depth[['Cohort','Patient ID','Sample ID','Median coverage']]
depth = depth[~depth['Sample ID'].str.contains('gDNA')]
depth = depth[~depth['Sample ID'].str.contains('_NL')]
depth = depth[~depth['Sample ID'].str.contains('cfDNA')]

tc = tc.merge(depth, on = ['Cohort','Sample ID','Patient ID'], how = 'left')
tc = tc[['Cohort','Patient ID','Sample ID','mut_TC','snp_TC','Median coverage']]

# =============================================================================
# Label low depth
# =============================================================================

tc['Color'] = 'k'
tc.loc[tc['Median coverage'] < 344.5, 'Color'] = 'orange'

# =============================================================================
# Scatter
# =============================================================================

fig,ax = plt.subplots()

lowdepth = tc[tc['Median coverage'] < 344.5]
highdepth = tc[tc['Median coverage'] > 344.5]

ax.scatter(lowdepth['mut_TC'], lowdepth['snp_TC'], color = lowdepth['Color'], s = 8, alpha = 0.45, clip_on = False, label = 'Depth < median')
ax.scatter(highdepth['mut_TC'], highdepth['snp_TC'], color = highdepth['Color'], s = 8, alpha = 0.45, clip_on = False, label = 'Depth > median')
ax.plot([0,1],[0,1], color = 'grey', linestyle = 'dashed')
ax.set_xlabel('mut_TC')
ax.set_ylabel('snp_TC')
ax.set_xlim(0,1)
ax.set_ylim(0,1)

tc = tc[tc['snp_TC'] != 0]
tc = tc[tc['mut_TC'] != 0]

lowdepth = tc[tc['Median coverage'] < 344.5]
highdepth = tc[tc['Median coverage'] > 344.5]

linregress_all = stats.linregress(tc['mut_TC'], tc['snp_TC'])
slope_all = str(round(linregress_all[0],2))
r_all = str(round(linregress_all[2],2))

linregress_high = stats.linregress(highdepth['mut_TC'], highdepth['snp_TC'])
slope_high = str(round(linregress_high[0],2))
r_high = str(round(linregress_high[2],2))


linregress_low = stats.linregress(lowdepth['mut_TC'], lowdepth['snp_TC'])
slope_low = str(round(linregress_low[0],2))
r_low = str(round(linregress_low[2],2))

ax.text(x = 0.02, y = 0.95, s = 'r_all: %s, r_high: %s, r_low: %s' % (r_all, r_high, r_low))
ax.legend()

plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/tumor_fraction_calcuation/TC_SNP_scatter_byDepth.pdf')


