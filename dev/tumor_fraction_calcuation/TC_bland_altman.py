# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:02:16 2021

@author: amurtha
"""

# =============================================================================
# Bland altman plot tumor content by category
# =============================================================================

import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# constants
# =============================================================================

color_dict = {1:'red',2:'blue',3:'red'}

# =============================================================================
# Import data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc['mut_TC'] = tc['mut_TC'].astype(float)
tc['snp_TC'] = tc['snp_TC'].astype(float)
tc['Group'] = tc['Group'].astype(int)

# =============================================================================
# Keep samples with an estimate from both
# =============================================================================

# tc = tc[(tc['mut_TC'] > 0)&(tc['snp_TC'] > 0)]

# =============================================================================
# Create columns
# =============================================================================

tc['Difference'] = tc['mut_TC'] - tc['snp_TC']
tc['Mean'] = (tc['mut_TC'] + tc['snp_TC']) / 2

# =============================================================================
# Set up color
# =============================================================================

tc['color'] = tc['Group'].replace(color_dict)
# tc = tc[tc['color'] == 'blue']

# =============================================================================
# Plot
# =============================================================================

fig,[ax1,ax2,ax3] = plt.subplots(nrows = 3, figsize = (7,7))

mean = tc['Difference'].mean()
std = tc['Difference'].std()

ax1.scatter(tc['Mean'],tc['Difference'], c = tc['color'], lw = 0, alpha = 0.5, label = 'Sample')
ax1.plot([0,1],[mean,mean], label = 'Bias')
ax1.plot([0,1],[mean+1.96*std,mean+1.96*std], ls = '--', c = 'orange', label = 'Limit of agreement')
ax1.plot([0,1],[mean-1.96*std,mean-1.96*std], ls = '--', c = 'orange')

ax1.set_xlim(0,1)
ax1.set_ylim(-1,1)
ax1.set_xlabel('Mean')
# ax1.set_ylabel('Difference between methods (Mut - SNP)')

ax1.legend(fontsize = 8)

ax2.scatter(tc['snp_TC'],tc['Difference'], c = tc['color'], lw = 0, alpha = 0.5, label = 'Sample')
ax2.plot([0,1],[mean,mean], label = 'Bias')
ax2.plot([0,1],[mean+1.96*std,mean+1.96*std], ls = '--', c = 'orange', label = 'Limit of agreement')
ax2.plot([0,1],[mean-1.96*std,mean-1.96*std], ls = '--', c = 'orange')

ax2.set_xlim(0,1)
ax2.set_ylim(-1,1)
ax2.set_xlabel('SNP TC')
ax2.set_ylabel('Difference between methods (Mut - SNP)')

ax3.scatter(tc['mut_TC'],tc['Difference'], c = tc['color'], lw = 0, alpha = 0.5, label = 'Sample')
ax3.plot([0,1],[mean,mean], label = 'Bias')
ax3.plot([0,1],[mean+1.96*std,mean+1.96*std], ls = '--', c = 'orange', label = 'Limit of agreement')
ax3.plot([0,1],[mean-1.96*std,mean-1.96*std], ls = '--', c = 'orange')

ax3.set_xlim(0,1)
ax3.set_ylim(-1,1)
ax3.set_xlabel('Mutation TC')
# ax3.set_ylabel('Difference between methods (Mut - SNP)')


plt.tight_layout()

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/tumor_fraction_calcuation/TC_bland_altman.pdf')

