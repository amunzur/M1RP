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
# Import data
# =============================================================================

neg_samples = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/abienza_panel_cna_controls.txt', header = None, names = ['Sample ID'])

cn = pd.read_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/Abi-enza_logratios_melted.tsv', sep = '\t',)


# =============================================================================
# Keep only confirmed negative samples
# =============================================================================

cn = cn[cn['Sample ID'].isin(neg_samples['Sample ID'])]

# =============================================================================
# Calculate distribution of copy number negative samples by gene
# =============================================================================

mean = cn.groupby('GENE').mean()['Log_ratio']
std = cn.groupby('GENE').std()['Log_ratio']
sem = cn.groupby('GENE').sem()['Log_ratio']
ci_diff = sem.apply(lambda x: x* stats.t.ppf((1+0.95)/2., len(cn['Sample ID'].unique())-1))

gene_noise = pd.DataFrame({'lr_mean':mean,'lr_std':std, 'lr_ci_low':mean-ci_diff, 'lr_ci_high':mean+ci_diff})

# =============================================================================
# Create graph
# =============================================================================

fig,ax = plt.subplots(figsize = (15,5))

m = 'o'
c = '#E6E7E8'
z = 100

ax.scatter(cn['GENE'],cn['Log_ratio'], marker = m, facecolor=c, edgecolor='k', s=50, alpha = 0.45, zorder = z)

for i, (index, row) in enumerate(gene_noise.reset_index().iterrows()):
    ax.plot([row['GENE'],row['GENE']], [row['lr_ci_low'],row['lr_ci_high']], lw = 1, zorder = 100, marker = '_', color = 'k')

ax.set_xticks(np.arange(0,len(cn['GENE'].unique()),1))
ax.set_xticklabels(cn['GENE'].unique(), rotation = 90)
    
ax.set_ylim(-1.5, 1.5)
ax.set_ylabel('Log ratio')
ax.set_xlim(-.7, len(cn['GENE'].unique())-0.3)
ax.grid(b=True, which='both', axis='y', color='0.65', linewidth=0.5, linestyle='dotted', zorder=0)
ax.grid(b=True, which='both', axis='x', color='0.65', linewidth=0.5, linestyle='dotted', zorder=0)
plt.tight_layout()


# =============================================================================
# Save figure and stats
# =============================================================================

gene_noise.to_csv('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/normal_samples_cn/abienza_normal_samples_cn.tsv', sep = '\t')
plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/normal_samples_cn/normal_samples_cn.png')
plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/CANdy/Abi_enza/normal_samples_cn/normal_samples_cn.pdf')