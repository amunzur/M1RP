# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 10:05:45 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as stats
import numpy as np
import sys


# =============================================================================
# Import data
# =============================================================================

if len(sys.argv) > 1:
    cohort = sys.argv[1]
else:
    cohort = 'M1RP'

mut_tc = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/tumor_fraction/%s_mut_tc.tsv' % cohort, sep = '\t')
snp_tc = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/tumor_fraction/%s_snp_tc.tsv' % cohort, sep = '\t')

snp_tc = snp_tc.merge(mut_tc, on = ['Cohort','Patient ID','Sample ID'])
snp_tc['Non-truncal flag'] = snp_tc['Non-truncal flag'].fillna(False)
snp_tc['snp_TC'] = snp_tc['snp_TC'].fillna(0)

snp_tc = snp_tc[['Cohort','Patient ID','Sample ID','mut_TC','Chromosome','Position', 
                 'Gene','Effect','Variant allele frequency','Read depth at variant position',
                 'Gene Log Ratio','Non-truncal flag','snp_TC', 'median LOH VAF','Gene count',
                 'LOH gene count','Copy neutral gene count']]

# =============================================================================
# Create final NGS tc
# =============================================================================

snp_tc['Final tNGS_TC'] = snp_tc['mut_TC']
snp_tc['Group'] = 1
snp_tc.loc[(snp_tc['Non-truncal flag'] == True), 'Final tNGS_TC'] = snp_tc['snp_TC']
snp_tc.loc[(snp_tc['Non-truncal flag'] == True), 'Group'] = 2

snp_tc.loc[(snp_tc['Variant allele frequency'] < 0.08)&(snp_tc['mut_TC'] > 0), 'Final tNGS_TC'] = snp_tc['snp_TC']
snp_tc.loc[(snp_tc['Variant allele frequency'] < 0.08)&(snp_tc['mut_TC'] > 0), 'Group'] = 2

snp_tc.loc[(snp_tc['mut_TC'] == 0)&(snp_tc['snp_TC'] > 0), 'Final tNGS_TC'] = snp_tc['snp_TC']
snp_tc.loc[(snp_tc['mut_TC'] == 0)&(snp_tc['snp_TC'] > 0), 'Group'] = 2


snp_tc.loc[(snp_tc['Group'] == 2)&(snp_tc['snp_TC'] == 0), 'Final tNGS_TC'] = snp_tc['mut_TC']
snp_tc.loc[(snp_tc['Group'] == 2)&(snp_tc['snp_TC'] == 0), 'Group'] = 3

snp_tc.loc[(snp_tc['mut_TC'] == 0)&(snp_tc['snp_TC'] == 0), 'Group'] = 4

# =============================================================================
# Create plots
# =============================================================================

fig = plt.figure(figsize = (4,3.5))
gs = gridspec.GridSpec(ncols=5, nrows=5, figure=fig, width_ratios = [0.2,1,0.4,0.2,1], height_ratios = [1,0.2,0.45,1,0.2])

ax1 = fig.add_subplot(gs[0,1])
ax1_ykde = fig.add_subplot(gs[0,0], sharey = ax1)
ax1_xkde = fig.add_subplot(gs[1,1], sharex = ax1)

ax1.set_zorder(10)
ax1_ykde.set_zorder(1)
ax1_xkde.set_zorder(1)

ax2 = fig.add_subplot(gs[0,4])
ax2_ykde = fig.add_subplot(gs[0,3], sharey = ax2)
ax2_xkde = fig.add_subplot(gs[1,4], sharex = ax2)

ax2.set_zorder(10)
ax2_ykde.set_zorder(1)
ax2_xkde.set_zorder(1)

ax3 = fig.add_subplot(gs[3,1])
'''
ax3_ykde = fig.add_subplot(gs[3,0])
ax3_xkde = fig.add_subplot(gs[4,1])
'''
ax4 = fig.add_subplot(gs[3, 4])


xpad = 13
pad = 16

# =============================================================================
# Plot ax1. Group1 snp_tc vs mut_tc
# =============================================================================

group1 = snp_tc[snp_tc['Group'] == 1]
ax1.scatter(group1['mut_TC'], group1['snp_TC'], s = 6, c = 'k', lw = 0, clip_on = False, alpha = 0.8)
ax1.plot([0,1],[0,1], c = 'gray', ls = 'dashed', zorder = 0, lw = 0.5)
ax1.set_ylim(0,1)
ax1.set_xlim(0,1)
ax1.set_xlabel('mut TC (%)', fontsize = 8)
ax1.set_ylabel('SNP TC (%)', fontsize = 8)
ax1.set_title('Truncal mutation used', fontsize = 8)
# group1 = group1[group1['snp_TC'] > 0]
r = stats.linregress(group1['mut_TC'],group1['snp_TC'])[2]
ax1.text(.05,0.9,'r = %s' % str(round(r,2)), ha ='left', fontsize = 8)

ax1.tick_params(labelsize = 8, pad = pad, width = 0.5)
ax1.tick_params(axis = "x", pad = xpad)

ax1.set_xticks(np.arange(0,1.1,0.25))
ax1.set_yticks(np.arange(0,1.1,0.25))

ax1.set_yticklabels(np.arange(0,110,25))
ax1.set_xticklabels(np.arange(0,110,25))

y_kde = stats.gaussian_kde(group1['snp_TC'])
x_kde = stats.gaussian_kde(group1['mut_TC'])

xs = np.linspace(0,1,100)
x_kde.covariance_factor = lambda : .25
x_kde._compute_covariance()

ax1_xkde.plot(xs,x_kde(xs), lw = 0.5)
ax1_xkde.fill_between(xs, 0, x_kde(xs))
ax1.set_ylim(bottom = 0)
ax1_xkde.invert_yaxis()

ax1_xkde.spines['left'].set_visible(False)
ax1_xkde.spines['bottom'].set_visible(False)
ax1_xkde.tick_params(left = False, bottom = False, labelleft=False, labelbottom=False)

ys = np.linspace(0,1,100)
y_kde.covariance_factor = lambda : .25
y_kde._compute_covariance()

ax1_ykde.plot(y_kde(ys), ys)
ax1_ykde.fill_betweenx(ys, 0, y_kde(ys))
ax1_ykde.set_xlim(left = 0)
ax1_ykde.invert_xaxis()

ax1_ykde.spines['left'].set_visible(False)
ax1_ykde.spines['bottom'].set_visible(False)
ax1_ykde.tick_params(left = False, bottom = False, labelleft=False, labelbottom=False)

# =============================================================================
# Plot ax2. Group2 snp_tc vs mut_tc
# =============================================================================

group2 = snp_tc[snp_tc['Group'] == 2]
ax2.scatter(group2['mut_TC'], group2['snp_TC'], s = 6, c = 'k', lw = 0, clip_on = False, alpha = 0.8)
ax2.plot([0,1],[0,1], c = 'gray', ls = 'dashed', zorder = 0, lw = 0.5)
ax2.set_ylim(0,1)
ax2.set_xlim(0,1)
ax2.set_xlabel('mut TC (%)', fontsize = 8)
ax2.set_title('No truncal mutations, SNPs used', fontsize = 8)
# group2 = group2[group2['mut_TC'] > 0]
r = stats.linregress(group2['mut_TC'],group2['snp_TC'])[2]
ax2.text(1,0.1,'r = %s' % str(round(r,2)), ha ='right', fontsize = 8)

ax2.tick_params(labelsize = 8, pad = pad, width = 0.5)
ax2.tick_params(axis = "x", pad = xpad)

ax2.set_xticks(np.arange(0,1.1,0.25))
ax2.set_yticks(np.arange(0,1.1,0.25))

ax2.set_yticklabels(np.arange(0,110,25))
ax2.set_xticklabels(np.arange(0,110,25))

y_kde = stats.gaussian_kde(group2['snp_TC'])
x_kde = stats.gaussian_kde(group2['mut_TC'])

xs = np.linspace(0,1,100)
x_kde.covariance_factor = lambda : .25
x_kde._compute_covariance()

ax2_xkde.plot(xs,x_kde(xs), lw = 0.5)
ax2_xkde.fill_between(xs, 0, x_kde(xs))
ax2.set_ylim(bottom = 0)
ax2_xkde.invert_yaxis()

ax2_xkde.spines['left'].set_visible(False)
ax2_xkde.spines['bottom'].set_visible(False)
ax2_xkde.tick_params(left = False, bottom = False, labelleft=False, labelbottom=False)

ys = np.linspace(0,1,100)
y_kde.covariance_factor = lambda : .25
y_kde._compute_covariance()

ax2_ykde.plot(y_kde(ys), ys)
ax2_ykde.fill_betweenx(ys, 0, y_kde(ys))
ax2_ykde.set_xlim(left = 0)
ax2_ykde.invert_xaxis()

ax2_ykde.spines['left'].set_visible(False)
ax2_ykde.spines['bottom'].set_visible(False)
ax2_ykde.tick_params(left = False, bottom = False, labelleft=False, labelbottom=False)

# =============================================================================
# Plot ax3. Group3 snp_tc vs mut_tc
# =============================================================================

group3 = snp_tc[snp_tc['Group'] == 3]
ax3.scatter(group3['mut_TC'], group3['snp_TC'], s = 6, c = 'k', lw = 0)
ax3.plot([0,1],[0,1], c = 'gray', ls = 'dashed', zorder = 0, lw = 0.5)
ax3.set_ylim(-0.02,1)
ax3.set_xlim(-0.02,1)
ax3.set_xlabel('mut TC', fontsize = 8)
ax3.set_ylabel('SNP TC', fontsize = 8)
ax3.set_title('Non-trunal mutation used, No SNP TC', fontsize = 8)
ax3.tick_params(labelsize = 8, width = 0.5)
ax3.set_xticks(np.arange(0,1.1,0.25))
ax3.set_yticks(np.arange(0,1.1,0.25))

ax3.set_xticklabels(np.arange(0,110,25))
ax3.set_yticklabels(np.arange(0,110,25))

'''y_kde = stats.gaussian_kde(group3['snp_TC'])
x_kde = stats.gaussian_kde(group3['mut_TC'])

xs = np.linspace(0,1,100)
x_kde.covariance_factor = lambda : .25
x_kde._compute_covariance()

ax3_xkde.plot(xs,x_kde(xs), lw = 0.5)
ax3_xkde.fill_between(xs, 0, x_kde(xs))
ax3_xkde.invert_yaxis()

ax3_xkde.spines['left'].set_visible(False)
ax3_xkde.spines['bottom'].set_visible(False)
ax3_xkde.tick_params(left = False, bottom = False, labelleft=False, labelbottom=False)

ys = np.linspace(0,1,100)
y_kde.covariance_factor = lambda : .25
y_kde._compute_covariance()

ax3_ykde.plot(y_kde(ys), ys)
ax3_ykde.fill_betweenx(ys, 0, y_kde(ys))
ax3_ykde.invert_xaxis()

ax3_ykde.spines['left'].set_visible(False)
ax3_ykde.spines['bottom'].set_visible(False)
ax3_ykde.tick_params(left = False, bottom = False, labelleft=False, labelbottom=False)'''

# =============================================================================
# Plot ax4. Swarm plots of all groups
# =============================================================================

snp_tc['x'] = snp_tc['Group'].apply(lambda x: x + np.random.rand()/5-1/10)

snp_tc_plot = snp_tc[snp_tc['Group'] != 4].copy()

ax4.scatter(snp_tc_plot['x'], snp_tc_plot['Final tNGS_TC'], s = 6, c = 'k', alpha = 0.3, lw = 0)
ax4.boxplot([group1['Final tNGS_TC'],group2['Final tNGS_TC'],group3['Final tNGS_TC']], sym = '',)
ax4.set_ylabel('Fianl tNGS TC', fontsize = 8)
ax4.set_xticks([1,2,3])
ax4.set_yticks(np.arange(0,1.1,0.25))
ax4.set_yticklabels(np.arange(0,110,25))
ax4.set_xticklabels(['Group1\nn=%i' % len(group1),'Group2\nn=%i' % len(group2),'Group3\nn=%i' % len(group3)])
ax4.tick_params(labelsize = 8, width = 0.5)

plt.tight_layout()
gs.update(hspace = 0, wspace = 0, bottom = 0.05, left = 0.10, right = 0.94)

plt.savefig('G:/Andy Murtha/Ghent/M1RP/prod/tumor_fraction_calcuation/M1RP_snpTC_mutTC_byGroup.pdf')
plt.savefig('G:/Andy Murtha/Ghent/M1RP/prod/tumor_fraction_calcuation/M1RP_snpTC_mutTC_byGroup.pdf')

snp_tc = snp_tc[['Cohort','Patient ID','Sample ID','mut_TC','Chromosome','Position', 'Gene','Effect','Variant allele frequency','Read depth at variant position','Gene Log Ratio','Non-truncal flag','Group','snp_TC', 'median LOH VAF','Gene count','LOH gene count','Copy neutral gene count','Final tNGS_TC']]


snp_tc.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/tumor_fraction/%s_tumor_fraction_final.tsv' % cohort, sep = '\t', index = None)

