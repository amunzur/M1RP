# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:55:18 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import random

'''
Plot: TC from 20-100 rolling standard deviation of SNP deviation in genes with 
> 3 snps and called hetz. loss
'''

# =============================================================================
# Import data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/M1RP_FFPE_cna.tsv', sep = '\t')
snps = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/SNPs/M1RP_snps_melted.tsv', sep = '\t')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc['Group'] = tc['Group'].astype(int)

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

# =============================================================================
# Keep only tissue
# =============================================================================

samples = samples[samples['Cohort'] == 'M1RP']
samples = samples[samples['Sample Category'].isin(['RP', 'MLN', 'PB', 'MB', 'MT'])]

tc = tc[tc['Sample ID'].isin(samples['Sample ID'].tolist())]

# =============================================================================
# Keep TC > 20. Remove Samples where tumor content was predicted by SNP deviation
# =============================================================================

tc = tc[tc['Final tNGS_TC'] > 0.2]
tc = tc[tc['Group'] != 2]

# =============================================================================
# Keep SNPs in called Hetz genes
# =============================================================================

cn = cn[cn['Copy_num'] < 0]
snps = snps.merge(cn[['Sample ID','GENE']], on = ['Sample ID','GENE'])

# =============================================================================
# Keep only genes with >=3 SNPs per gene. Groupby and median
# =============================================================================

snps['Count'] = 1
snps_count = snps[['Sample ID','GENE','Count']].groupby(['Sample ID','GENE']).count()
snps_count = snps_count[snps_count['Count'] >= 3].reset_index()

snps = snps.drop('Count', axis = 1)
snps = snps.merge(snps_count, on = ['Sample ID', 'GENE'])

snps = snps.groupby(['Sample ID','GENE']).median().reset_index()
snps['VAF'] = (snps['VAF']-0.5).abs()
snps = snps.merge(tc[['Sample ID','Final tNGS_TC']], on = 'Sample ID')

# =============================================================================
# Create plot
# =============================================================================

fig,[ax1,ax2] = plt.subplots(ncols = 2)

# =============================================================================
# Create rolling average
# =============================================================================

xs = []
ys = []

for x in np.arange(0.2, 1.0, 0.05):
    tmp = snps[(snps['Final tNGS_TC'] < x+0.05)&(snps['Final tNGS_TC'] >= x-0.05)].copy()
    ax1.scatter(x,tmp['VAF'].mean(), c = 'k')
    xs.append(x)
    ys.append(tmp['VAF'].mean())
    
lin = stats.linregress(xs,ys)

ax1.plot(np.arange(0.2,1,0.01),(1/(2-np.arange(0.2,1,0.01))-0.5), label = 'Theoretical deviation')

ax1.set_xlabel('Tumor content (rolling average +/-0.05)')
ax1.set_ylabel('Mean absolute SNP deviation')
ax1.legend(fontsize = 8)

ax1.text(0.2,0.45,'r = %.3f' % lin[2])

# =============================================================================
# Plot boxplots and swarms
# =============================================================================

med_tc = snps.drop_duplicates('Sample ID')['Final tNGS_TC'].median()

above = snps[snps['Final tNGS_TC'] >= med_tc].copy()
below = snps[snps['Final tNGS_TC'] <  med_tc].copy()

above['x'] = 2
above['x'] = above['x'].apply(lambda x: x + random.uniform(-0.4, 0.4))

below['x'] = 1
below['x'] = below['x'].apply(lambda x: x + random.uniform(-0.4, 0.4))

# ax2.boxplot([below['Log_ratio'].tolist(),above['Log_ratio'].tolist()], showfliers = False)
ax2.scatter(above['x'],above['VAF'], s = 2, lw = 0, c = 'k', alpha = 0.5)
ax2.scatter(below['x'],below['VAF'], s = 2, lw = 0, c = 'k', alpha = 0.5)

ax2.text(0.6, 0.45, 'Mann-Whitney U (1 sided)\np = %.3f\nTC median = %.3f' % (stats.mannwhitneyu(above['VAF'], below['VAF'], alternative='greater')[1], med_tc))

ax2.set_xticks([1,2])
ax2.set_xticklabels(['TC < median', 'TC >= median'])
ax2.set_ylim(0,0.5)
ax2.set_ylabel('Mean absolute SNP deviation')

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/summary/SNP_vs_TC.pdf')