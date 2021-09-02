# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:02:16 2021

@author: amurtha
"""

# =============================================================================
# Bland altman plot tumor content by category
# =============================================================================

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import gridspec

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

tc = tc[tc['Group'] != 4]

# =============================================================================
# Create columns
# =============================================================================

tc['Difference'] = tc['mut_TC'] - tc['snp_TC']
tc['Mean'] = (tc['mut_TC'] + tc['snp_TC']) / 2

# =============================================================================
# Set up color
# =============================================================================

tc['color'] = tc['Group'].replace(color_dict)
tc.loc[tc['snp_TC'] == 0, 'color'] = 'lightgrey' 
tc.loc[tc['mut_TC'] == 0, 'color'] = 'lightgrey'

# =============================================================================
# Plot
# =============================================================================

fig = plt.figure(figsize = (4,3.5))
gs = gridspec.GridSpec(ncols = 2, nrows = 2, figure = fig)

ax1 = fig.add_subplot(gs[0,0:2])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[1,1])

mean = tc[(tc['snp_TC'] != 0)&(tc['mut_TC'] != 0)]['Difference'].mean()
std = tc[(tc['snp_TC'] != 0)&(tc['mut_TC'] != 0)]['Difference'].std()

ax1.scatter(tc['Mean'],tc['Difference'], c = tc['color'], lw = 0, alpha = 0.5, s = 12)
ax1.plot([0,1],[mean,mean], label = 'Bias', lw = 0.8, c = 'k')
ax1.plot([0,1],[mean+1.96*std,mean+1.96*std], ls = '--', c = 'orange', label = 'Limit of agreement', lw = 0.8)
ax1.plot([0,1],[mean-1.96*std,mean-1.96*std], ls = '--', c = 'orange', lw = 0.8)

ax1.plot([0,0.5],[0,1], color = 'black', ls = 'dashed', lw = 0.5, zorder = 0)
ax1.plot([0,0.5],[0,-1], color = 'black', ls = 'dashed', lw = 0.5, zorder = 0)

ax1.set_xlim(0,1)
ax1.set_ylim(-1,1)
ax1.set_xlabel('Mean', fontsize = 6)
ax1.set_ylabel('Difference b/w methods\n(Mut - SNP)', fontsize = 6)

handles = [Line2D([0],[0], lw = 0, color = 'r', marker = 'o', markeredgewidth=0, markersize = 4),
           Line2D([0],[0], lw = 0, color = 'b', marker = 'o', markeredgewidth=0, markersize = 4),
           Line2D([0],[0], lw = 0.8, color = 'k', marker = None , markeredgewidth=0),
           Line2D([0],[0], lw = 0.8, color = 'orange', marker = None , markeredgewidth=0),]
labels = ['SNP TC used','Mut. TC used','Bias','Limit of agreement']

ax1.legend(handles, labels, fontsize = 6)
ax1.tick_params(labelsize = 6)

# =============================================================================
# Plot differences histogram
# =============================================================================

dist = tc[(tc['snp_TC'] != 0)&(tc['mut_TC'] != 0)]['Difference']

ax2.hist(dist, bins = 15)
mean = dist.mean(); std = dist.std(); variance = np.square(std)
x = np.arange(-1,1,.001)

ax2.plot(x, stats.norm.pdf(x, mean,std)*30)
ax2.set_ylabel('Frequency', fontsize = 6)
ax2.set_xlabel('Difference (Mut - SNP)', fontsize = 6)

ax2.tick_params(width=0.5, labelsize = 6)

sw = stats.shapiro(dist)
ax2.text(1.0, 70, "Shapiro-Wilk\np = %.5f" % sw[1], fontsize = 6, ha = 'right')

# =============================================================================
# Ax3 QQ plot for normal scores
# =============================================================================

data, r = stats.probplot(dist, dist='norm', fit = True)
ax3.scatter(data[0],data[1], s = 10, lw = 0, alpha = 0.6, zorder = 10)
ax3.plot([-3,3], [r[0]*-3+r[1],r[0]*3+r[1]], ls = 'dashed', lw = 0.5, color = 'k', zorder = 100)
r_s2 = r[2]**2

ax3.text(-2.5, 0.4, "$R^{2}$ = %0.3f" % r_s2, fontsize = 6, ha = 'left')

ax3.tick_params(labelsize = 6, width = 0.5)
ax3.set_ylabel('Ordered Values', fontsize = 6)
ax3.set_xlabel('Theoretical quantiles', fontsize = 6)



# ax2.set_ylabel('Difference between methods (Mut - SNP)')
'''
ax3.scatter(tc['mut_TC'],tc['Difference'], c = tc['color'], lw = 0, alpha = 0.5, label = 'Sample')
ax3.plot([0,1],[mean,mean], label = 'Bias')
ax3.plot([0,1],[mean+1.96*std,mean+1.96*std], ls = '--', c = 'orange', label = 'Limit of agreement')
ax3.plot([0,1],[mean-1.96*std,mean-1.96*std], ls = '--', c = 'orange')

ax3.set_xlim(0,1)
ax3.set_ylim(-1,1)
ax3.set_xlabel('Mutation TC')
# ax3.set_ylabel('Difference between methods (Mut - SNP)')
'''

plt.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/TC_bland_altman.pdf')
fig.savefig('G:/Andy Murtha/Ghent/M1RP/prod/tumor_fraction_calcuation/TC_bland_altman.pdf')
