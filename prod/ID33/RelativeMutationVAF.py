# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 12:25:02 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import string
import numpy as np
from matplotlib.lines import Line2D

color_dict = {('TP53','Splice'):'blue',
              ('RNF43','Splice'):'red',
              ('PTEN','Non-frameshift'):'green',
              ('PTEN','Frameshift'):'yellow'}
sample_type_dict = {'PB':0,'RP':1,'MLN':2,'MB':2,'MT':2,'cfDNA':3}


# =============================================================================
# Import mutation TC and CN data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

# =============================================================================
# Keep ID33
# =============================================================================

muts = muts[muts['Patient ID'] == 'ID33']
tc = tc[tc['Patient ID'] == 'ID33']

# =============================================================================
# Grab sample with mutation in these genes. Label color of the mutations
# =============================================================================

muts = muts[muts['GENE'].isin(['BRCA2','TP53','RNF43','PTEN'])]
muts['Color'] = muts.apply(lambda x: color_dict.get((x['GENE'],x['EFFECT'].split(' ')[0])), axis = 1)

# =============================================================================
# Merge tc. ORder by sample type then TC
# =============================================================================
muts['Sample category'] = muts['Sample ID'].str.split('_').str[2].str.strip(string.digits)
muts['Sample category'] = muts['Sample category'].map(sample_type_dict)

# =============================================================================
# Create relative VAF
# =============================================================================

brca2_vaf = muts.copy()
brca2_vaf = brca2_vaf[brca2_vaf['GENE'] == 'BRCA2']
brca2_vaf = brca2_vaf[['Sample ID','Allele_frequency']]
brca2_vaf.columns = ['Sample ID','brca2_vaf']

muts = muts[muts['GENE'] != 'BRCA2']
muts = muts.merge(brca2_vaf, on = 'Sample ID')

muts['Relative VAF'] = muts['Allele_frequency'] - muts['brca2_vaf']

# =============================================================================
# Split by number of mutations present in the sample
# =============================================================================

grouped = muts.copy().groupby(['Sample ID']).count()
grouped = grouped[['GENE']].reset_index()
grouped.columns = ['Sample ID','Count']

muts = muts.merge(grouped, on ='Sample ID')

muts['Sample Category'] = muts['Sample ID'].str.split('_').str[2].str.strip(string.digits)
muts = muts.sort_values(['Sample Category','Count','Final tNGS_TC'], ascending = [False,True, False])

# =============================================================================
# Create figure
# =============================================================================

fig,ax = plt.subplots(figsize = (2.35,2.15))

# =============================================================================
# Plot 
# =============================================================================

ax.scatter(muts['Sample ID'],muts['Relative VAF'], c = muts['Color'], lw = 0, s = 12,zorder = 10)
ax.plot([-1,len(muts['Sample ID'].drop_duplicates())],[0,0], ls = 'dashed', lw = 0.8, color = 'lightgrey', zorder = 0)

ax.set_xticks(np.arange(0,len(muts['Sample ID'].drop_duplicates()),1))
ax.set_xticklabels(muts['Sample ID'].drop_duplicates().str.split('_').str[2], rotation = 90, fontsize = 8)

ax.set_xlim(-0.6, len(muts['Sample ID'].drop_duplicates())-0.4)

ax.set_ylim(-0.4,0.4)
ax.set_ylabel('Relative VAF', fontsize = 6, labelpad = 0)

fig.legend([Line2D([0],[0], marker = 'o', lw = 0, markeredgewidth=0, markerfacecolor=x, markersize = 4) for x in list(color_dict.values())], [x for x in [" ".join(y) for y in list(color_dict.keys())]], fontsize = 6, ncol = 2, columnspacing = 1.6, handletextpad = 0.5)

ax.tick_params(labelsize = 6, bottom = False, pad = 0, size = 3)


plt.tight_layout()
fig.subplots_adjust(top = 0.98, left = 0.15, right = 0.97, bottom = 0.17)

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Patient Specific/ID33/relative VAF.pdf', transparent=True)