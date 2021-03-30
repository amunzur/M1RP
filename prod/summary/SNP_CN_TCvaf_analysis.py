# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 11:18:29 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random

# =============================================================================
# Constants
# =============================================================================

cn_color = {-2:'#3F60AC',
            -1:'#9CC5E9',
            0: '#EFEFEF',
            1: '#F59496',
            2: '#EE2D24'}

# =============================================================================
# Import data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
snps = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/hetz_snps/M1RP_snps_melted.tsv', sep = '\t')
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
samples = samples[samples['Sample Category'].isin(['RP', 'MLN', 'PB', 'NL', 'MB', 'MT', 'UCC', 'NT'])]

cn = cn[cn['Sample ID'].isin(samples['Sample ID'].tolist())]
muts = muts[muts['Sample ID'].isin(samples['Sample ID'].tolist())]
snps = snps[snps['Sample ID'].isin(samples['Sample ID'].tolist())]

samples = samples[['Sample ID','Sample Category']].merge(tc[['Sample ID','Group','Final tNGS_TC']], on = 'Sample ID', how = 'left')
samples['Final tNGS_TC'] = samples['Final tNGS_TC'].fillna(0)
samples.loc[(samples['Final tNGS_TC'] == 0)&(~samples['Group'].isnull()), 'Group'] = 4
samples['Group'] = samples['Group'].fillna(5)

# =============================================================================
# Prepare patient order
# Order first by TC calculation type, then TC
# =============================================================================

sample_order = samples.copy().sort_values(['Group','Final tNGS_TC'], ascending = [True,False])['Sample ID'].tolist()

groups = samples[['Sample ID','Group']].groupby('Group').count()
groups.columns = ['Group count']

# This code will add a 5 row spacer between groups
s_order = pd.DataFrame({'Sample ID':sample_order, 'y': np.arange(0,len(sample_order),1)})
s_order = s_order.merge(samples[['Sample ID','Group']], on = 'Sample ID')
s_order['y'] = s_order['y'] + (s_order['Group']-1) * 5

s_order = dict(zip(s_order['Sample ID'].tolist(),s_order['y'].tolist()))

# =============================================================================
# Label y-coordinates for tumor fraction and mutations. 
# =============================================================================

samples['y'] = samples['Sample ID'].apply(lambda x: s_order.get(x))
muts['y'] = muts['Sample ID'].apply(lambda x: s_order.get(x))
cn['y'] = cn['Sample ID'].apply(lambda x: s_order.get(x))

# =============================================================================
# Limit cn analysis to certain genes. Add color column
# =============================================================================

bad_genes = ['CDKN1B','AKT1','ERCC2','AKT3','HSD3B1','TMPRSS2']

cn['Log_ratio'] = cn['Log_ratio'].abs()
cn['Color'] = cn['Copy_num'].apply(lambda x: cn_color.get(x))
cn = cn[~cn['GENE'].isin(bad_genes)]


# =============================================================================
# Group SNPs by gene, sample ID. Keep only genes with >3 SNPs that are CN or 1 copy loss?
# =============================================================================

cn_snp = cn.copy()
cn_snp = cn_snp[cn_snp['Copy_num'] <= 0]
cn_snp = cn_snp[['Sample ID','GENE']]

snps['Count'] = 1
snps_count = snps[['Sample ID','GENE','Count']].groupby(['Sample ID','GENE']).count()
snps_count = snps_count[snps_count['Count'] >= 3].reset_index()

snps = snps.drop('Count', axis = 1)
snps = snps.merge(snps_count, on = ['Sample ID', 'GENE'])
snps = snps.merge(cn_snp, on = ['Sample ID','GENE'])

snps = snps[['Sample ID','GENE','VAF']].groupby(['Sample ID','GENE']).median().reset_index()
snps['y'] = snps['Sample ID'].apply(lambda x: s_order.get(x))
snps['VAF'] = (snps['VAF'] - 0.5).abs() + 0.5

# =============================================================================
# Create axis
# =============================================================================

fig,[ax,g0,ax1,g1,ax2,g2,ax3] = plt.subplots(ncols = 7, sharey = True, figsize = (7.5,7), gridspec_kw={'width_ratios':[0.03,0.02,1,0.1, 1,0.1, 1]})

g0.set_visible(False)
g1.set_visible(False)
g2.set_visible(False)

# =============================================================================
# Plot tumor fraction and mutation
# =============================================================================

ax1.barh(samples['y'], samples['Final tNGS_TC'], zorder = 10, label = 'Tumor content', color = 'grey')
ax1.scatter(muts['Allele_frequency'],muts['y'], c = 'k', s = 7, alpha = 0.8, zorder = 100, label = 'Mutation VAF', clip_on = False, lw = 0)

ax1.set_ylim(-1,max(s_order.values()))
ax1.set_xlim(0,1)
ax1.set_xlabel('Tumor content (VAF)')
ax1.tick_params(left = False, labelleft = False)
ax1.legend(fontsize = 8)

# =============================================================================
# Plot absolute CN log-ratio
# =============================================================================

ax2.scatter(cn['Log_ratio'], cn['y'], s = 7, c = cn['Color'], alpha = 0.5, clip_on = False, lw = 0)

ax2.set_xlim(-0.03,3)
ax2.set_xlabel('Absolute copy-number\nlog-ratio')

# =============================================================================
# Plot average SNP allele frequency
# =============================================================================

ax3.scatter(snps['VAF'],snps['y'], s = 7, c = 'k', alpha = 0.25, clip_on = False, lw = 0)
# ax3.hexbin(snps['VAF'],snps['y'], gridsize = 100, bins='log', mincnt = 2)

ax3.set_xlim(0.495,1)
ax3.set_xticks(np.arange(0.5,1.01,0.1))
ax3.set_xlabel('Median SNP deviation\n(No LR >=0.3)')

# =============================================================================
# Plot Groups onto Y axis
# =============================================================================

bottom = 0
ys = []
for index,row in groups.iterrows():
    ax.bar(0,row['Group count']-0.5, bottom = bottom - 0.5)
    ys.append(bottom+(row['Group count'] / 2))
    bottom = bottom + row['Group count'] + 5
    
    
ax.set_yticks(ys)
ax2.set_yticks(ys)
ax3.set_yticks(ys)
ax.set_yticklabels(['Mutation-based TC','SNP-based TC','Non-truncal\nmut. TC','TC-neg.','NT'], rotation = 90, va = 'center')
ax1.set_ylim(-1,max(s_order.values()))
    
ax.tick_params(left = False,bottom = False, labelbottom = False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)    

    
# =============================================================================
# Adjust and save    
# =============================================================================

fig.tight_layout()
fig.subplots_adjust(wspace = 0)

fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/summary/SNP_CN_TCvaf_analysis.pdf')
fig.savefig('G:/Andy Murtha/Ghent/M1RP/dev/summary/SNP_CN_TCvaf_analysis.png', dpi = 500)