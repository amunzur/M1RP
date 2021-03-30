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

snp_tc = snp_tc[['Cohort','Patient ID','Sample ID','mut_TC','Chromosome','Position', 'Gene','Effect','Variant allele frequency','Read depth at variant position','Gene Log Ratio','Non-truncal flag','snp_TC', 'median LOH VAF','Gene count','LOH gene count','Copy neutral gene count']]

# =============================================================================
# Create final NGS tc
# =============================================================================

snp_tc['Final tNGS_TC'] = snp_tc['mut_TC']
snp_tc['Group'] = 1
snp_tc.loc[(snp_tc['Non-truncal flag'] == False), 'Final tNGS_TC'] = snp_tc['snp_TC']
snp_tc.loc[(snp_tc['Non-truncal flag'] == False), 'Group'] = 2

snp_tc.loc[(snp_tc['Variant allele frequency'] < 0.08)&(snp_tc['mut_TC'] > 0), 'Final tNGS_TC'] = snp_tc['snp_TC']
snp_tc.loc[(snp_tc['Variant allele frequency'] < 0.08)&(snp_tc['mut_TC'] > 0), 'Group'] = 2


snp_tc.loc[(snp_tc['Group'] == 2)&(snp_tc['snp_TC'] == 0), 'Final tNGS_TC'] = snp_tc['mut_TC']
snp_tc.loc[(snp_tc['Group'] == 2)&(snp_tc['snp_TC'] == 0), 'Group'] = 3

# =============================================================================
# Create plots
# =============================================================================

fig = plt.figure(constrained_layout=True)
gs = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

# =============================================================================
# Plot ax1. Group1 snp_tc vs mut_tc
# =============================================================================

group1 = snp_tc[snp_tc['Group'] == 1]
ax1.scatter(group1['mut_TC'], group1['snp_TC'], s = 3, c = 'k')
ax1.plot([0,1],[0,1], c = 'gray', ls = 'dashed', zorder = 0)
ax1.set_ylim(-0.02,1)
ax1.set_xlim(-0.02,1)
ax1.set_xlabel('mut TC')
ax1.set_ylabel('SNP TC')
ax1.set_title('Mutation "truncal", Mutation used', fontsize = 10)
# group1 = group1[group1['snp_TC'] > 0]
r = stats.linregress(group1['mut_TC'],group1['snp_TC'])[2]
ax1.text(.05,0.9,'r = %s' % str(round(r,2)), ha ='left')

# =============================================================================
# Plot ax2. Group2 snp_tc vs mut_tc
# =============================================================================

group2 = snp_tc[snp_tc['Group'] == 2]
ax2.scatter(group2['mut_TC'], group2['snp_TC'], s = 3, c = 'k')
ax2.plot([0,1],[0,1], c = 'gray', ls = 'dashed', zorder = 0)
ax2.set_ylim(-0.02,1)
ax2.set_xlim(-0.02,1)
ax2.set_xlabel('mut TC')
ax2.set_ylabel('SNP TC')
ax2.set_title('Mutation "nontruncal", SNP used', fontsize = 10)
# group2 = group2[group2['mut_TC'] > 0]
r = stats.linregress(group2['mut_TC'],group2['snp_TC'])[2]
ax2.text(1,0.1,'r = %s' % str(round(r,2)), ha ='right')

# =============================================================================
# Plot ax3. Group3 snp_tc vs mut_tc
# =============================================================================

group3 = snp_tc[snp_tc['Group'] == 3]
ax3.scatter(group3['mut_TC'], group3['snp_TC'], s = 3, c = 'k')
ax3.plot([0,1],[0,1], c = 'gray', ls = 'dashed', zorder = 0)
ax3.set_ylim(-0.02,1)
ax3.set_xlim(-0.02,1)
ax3.set_xlabel('mut TC')
ax3.set_ylabel('SNP TC')
ax3.set_title('Mutation "nontruncal", No SNP TC\nMutation used', fontsize = 10)

# =============================================================================
# Plot ax4. Swarm plots of all groups
# =============================================================================

snp_tc['x'] = snp_tc['Group'].apply(lambda x: x + np.random.rand()/5-1/10)

ax4.scatter(snp_tc['x'], snp_tc['Final tNGS_TC'], s = 3, c = 'k', alpha = 0.4)
ax4.boxplot([group1['Final tNGS_TC'],group2['Final tNGS_TC'],group3['Final tNGS_TC']], sym = '',)
ax4.set_ylabel('Fianl tNGS TC')
ax4.set_xticks([1,2,3])
ax4.set_xticklabels(['Group1','Group2','Group3'])

plt.tight_layout()
plt.savefig('G:/Andy Murtha/Ghent/M1RP/prod/tumor_fraction_calcuation/%s_snpTC_mutTC_byGroup.pdf' % cohort)

snp_tc = snp_tc[['Cohort','Patient ID','Sample ID','mut_TC','Chromosome','Position', 'Gene','Effect','Variant allele frequency','Read depth at variant position','Gene Log Ratio','Non-truncal flag','Group','snp_TC', 'median LOH VAF','Gene count','LOH gene count','Copy neutral gene count','Final tNGS_TC']]


snp_tc.to_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/tumor_fraction/%s_tumor_fraction_final.tsv' % cohort, sep = '\t', index = None)

