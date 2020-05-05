# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 13:20:03 2020

@author: amurtha
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats as stats

# =============================================================================
# Helpers
# =============================================================================


# =============================================================================
# Import data
# =============================================================================

samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/gene_cna_melted.tsv', sep = '\t')

# =============================================================================
# Keep FFPE samples, sample with low (< 10%) TC
# =============================================================================

sample_list_to_remove = [
    'M1RP_ID12_PB2',
    'M1RP_ID12_RP1',
    'M1RP_ID12_RP4',
    'M1RP_ID12_RP6',
    'M1RP_ID12_RP7',
    'M1RP_ID12_RP8',
    'M1RP_ID14_RP9']

samples['tNGS_TC'] = samples['tNGS_TC'].str.split('%').str[0].astype(np.float64) / 100.
samples['Sample type'] = samples['Sample ID'].str.split('_').str[2].str.extract('(\D+)')
samples = samples[samples['Sample type'].isin(['RP', 'MLN', 'PB'])]
samples = samples[samples['tNGS_TC'] < .10]
samples = samples[~samples['Sample ID'].isin(sample_list_to_remove)]

cn = cn[(cn['Sample ID'].isin(samples['Sample ID']))|(cn['Sample ID'].str.contains('NL'))]

# =============================================================================
# Calculate normal cell content
# =============================================================================

samples['tNGS NC'] = 1 - samples['tNGS_TC']

# =============================================================================
# Correct copy number for normal content
# =============================================================================

cn = cn.merge(samples[['Sample ID','tNGS NC']], on = 'Sample ID', how = 'left')
#Add 100% NC for NL samples
cn.loc[cn['Sample ID'].str.contains('NL'), 'tNGS NC'] = 1

cn['copies'] = -1
cn.loc[cn['CHROMOSOME'] == 'chrX','copies']=1*(1 + (2**cn['log_ratio']-1)/cn['tNGS NC'])
cn.loc[cn['CHROMOSOME'] != 'chrX','copies']=2*(1 + (2**cn['log_ratio']-1)/cn['tNGS NC'])
cn.loc[cn['copies'] < 0, 'copies'] = 0

# =============================================================================
# Calculate per gene stats to be plotted and saved
# =============================================================================

mean = cn.groupby('GENE').mean()['log_ratio']
std = cn.groupby('GENE').std()['log_ratio']
sem = cn.groupby('GENE').sem()['log_ratio']
ci_diff = sem.apply(lambda x: x* stats.t.ppf((1+0.95)/2., len(cn['Sample ID'].unique())-1))

gene_noise = pd.DataFrame({'lr_mean':mean,'lr_std':std, 'lr_ci_low':mean-ci_diff, 'lr_ci_high':mean+ci_diff})

mean = cn.groupby('GENE').mean()['copies']
std = cn.groupby('GENE').std()['copies']
sem = cn.groupby('GENE').sem()['copies']
ci_diff = sem.apply(lambda x: x* stats.t.ppf((1+0.95)/2., len(cn['Sample ID'].unique())-1))

gene_noise['copy_mean'] = mean
gene_noise['copy_std'] = std
gene_noise['copy_ci_low'] = mean-ci_diff
gene_noise['copy_ci_high'] = mean+ci_diff

# =============================================================================
# Create graph
# =============================================================================

fig,ax = plt.subplots(figsize = (15,5))

m = 'o'
c = '#E6E7E8'
z = 100

ax.scatter(cn['GENE'],cn['copies'], marker = m, facecolor=c, edgecolor='k', s=50, alpha = 0.45, zorder = z)

for i, (index, row) in enumerate(gene_noise.reset_index().iterrows()):
    ax.plot([row['GENE'],row['GENE']], [row['copy_ci_low'],row['copy_ci_high']], lw = 1, zorder = 100, marker = '_', color = 'k')

ax.set_xticks(np.arange(0,len(cn['GENE'].unique()),1))
ax.set_xticklabels(cn['GENE'].unique(), rotation = 90)
    
ax.set_ylim(-.2, 5)
ax.set_ylabel('Number of copies')
ax.set_xlim(-.7, len(cn['GENE'].unique())-0.3)
ax.grid(b=True, which='both', axis='y', color='0.65', linewidth=0.5, linestyle='dotted', zorder=0)
ax.grid(b=True, which='both', axis='x', color='0.65', linewidth=0.5, linestyle='dotted', zorder=0)
plt.tight_layout()


# =============================================================================
# Save figure and stats
# =============================================================================

# gene_noise.to_csv('G:\\Andy Murtha\\Ghent\M1RP\\dev\\normal_samples_cn\\normal_samples_cn.tsv', sep = '\t')
# plt.savefig('G:\\Andy Murtha\\Ghent\M1RP\\dev\\Figures\\normal_samples_cn\\normal_samples_cn.pdf')