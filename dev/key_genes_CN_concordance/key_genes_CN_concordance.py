# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 10:40:24 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Constants
# =============================================================================

key_genes = ['TP53','PTEN','RB1','CHD1','CDKN2A','AR']

# =============================================================================
# import data
# =============================================================================

cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/M1RP_FFPE_cna.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100

# =============================================================================
# Keep only samples with > 20% tc, merge tc onto cn
# =============================================================================

tc = tc[tc['Final tNGS_TC'] >= 0.4]
cn = cn[cn['Sample ID'].isin(tc['Sample ID'].tolist())]
tc = tc[['Sample ID','Final tNGS_TC']]
cn = cn.merge(tc, on = 'Sample ID')

# =============================================================================
# Add color
# =============================================================================

cn['Color'] = 'grey'
cn.loc[(cn['Sample ID'].str.contains('MLN'))|(cn['Sample ID'].str.contains('MB'))|(cn['Sample ID'].str.contains('cfDNA')), 'Color'] = 'red'

# =============================================================================
# Calculate absolute CN
# =============================================================================

cn['ploidy'] = 2
cn.loc[cn['CHROMOSOME'] == 'chrX', 'ploidy'] = 1

cn['Copies'] = (cn['ploidy'] * (1 + (2**(cn['Log_ratio'])-1)/cn['Final tNGS_TC'])).clip(lower = 0)

# =============================================================================
# Plot figures
# =============================================================================

fig, axs = plt.subplots(nrows = len(key_genes), sharex = True, figsize = (7.5,10))

for (ax,gene) in zip(axs,key_genes):
    cn_gene = cn[cn['GENE'] == gene]
    ax.scatter(cn_gene['Patient ID'], cn_gene['Copies'], c = cn_gene['Color'], alpha = 0.5, s = 6, zorder = 100)
    ax.set_xlabel(gene)
    ax.set_ylim(0,4)
    ax.set_yticks(np.arange(0,5,1))
    ax.set_xticklabels(cn_gene['Patient ID'].unique(), rotation = 90)
    ax.grid(axis = 'y', linestyle = 'dashed', color = '#dedede', zorder = 0)
    ax.tick_params(bottom = False)
    
plt.tight_layout()
plt.savefig('G:/Andy Murtha/Ghent/M1RP/dev/key_genes_CN_concordance/key_genes_CN_concordance.pdf')