# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 19:42:56 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import string
import scipy.stats as stats
import numpy as np
from matplotlib.patches import Patch 

# =============================================================================
# Constants
# =============================================================================

s_type_dict = {'PB':'Primary',
               'RP':'Primary',
               'MLN':'Metastatic',
               'MB':'Metastatic',
               'MT':'Metastatic',
               'cfDNA':'cfDNA',}

mut_cat_dict = {'Missense':'Missense',
                "5'-UTR":'Silent',
                'Intronic':'Silent',
                'Intergenic':'Silent',
                'Frameshift':'Truncation',
                'Stopgain':'Truncation',
                'Synonymous':'Silent',
                "3'-UTR":'Silent',
                'Non-frameshift':'Non-frameshift',
                'Exonic':'Silent',
                'Splice':'Truncation'}

mut_color_dict = {'Missense':'green',
                  'Truncation':'yellow',
                  'Non-frameshift':'black',
                  'Silent':'grey'}

# =============================================================================
# Helpers
# =============================================================================

    
# =============================================================================
# Main: 
#   Import mutations and tumor purity data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
bet = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/betastasis/M1RP_betastasis_all.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])

tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

# =============================================================================
# Limit to TC+ samples
# =============================================================================

tc = tc[tc['Final tNGS_TC'] > 0.1]
tc = tc[tc['Patient ID'] != "ID8"]

tc['Sample type'] = tc['Sample ID'].apply(lambda x: s_type_dict.get(x.split('_')[2].rstrip(string.digits)))

# =============================================================================
# Merge TC onto muts
# =============================================================================

muts = muts.merge(tc[['Sample ID','Sample type','Final tNGS_TC']], how = 'left', on = 'Sample ID')
muts['aVAF'] = muts['Allele_frequency'] / muts['Final tNGS_TC']

# =============================================================================
# Split by primary and metastatic
# =============================================================================

muts = muts[muts['Sample type'].isin(['Primary','Metastatic'])]
muts = muts.set_index(['Patient ID','GENE','POSITION','EFFECT'])

p_muts = muts[muts['Sample type'] == 'Primary'].copy()
m_muts = muts[muts['Sample type'] == 'Metastatic'].copy()

# =============================================================================
# Find percent exclusive and shared
# =============================================================================

print("Total unique muts: %i" % (len(muts.index.drop_duplicates())))

p_ex_muts = p_muts[~p_muts.index.isin(m_muts.index)]
m_ex_muts = m_muts[~m_muts.index.isin(p_muts.index)]

shared_muts = p_muts[p_muts.index.isin(m_muts.index)]

print("Shared mutations: %i" % (len(shared_muts.index.drop_duplicates())))
print("P_ex mutations: %i" % (len(p_ex_muts.index.drop_duplicates())))
print("M_ex mutations: %i" % (len(m_ex_muts.index.drop_duplicates())))

# =============================================================================
# Relabel back onto muts
# =============================================================================

muts['mut_loc'] = 'p_ex'
muts.loc[muts.index.isin(m_ex_muts.index), 'mut_loc'] = 'm_ex'
muts.loc[muts.index.isin(shared_muts.index), 'mut_loc'] = 'shared'

# =============================================================================
# Bar plot of mutation count by mutation type
# =============================================================================

muts = muts.reset_index()
muts['mut_cat'] = muts['EFFECT'].apply(lambda x: mut_cat_dict.get(x.split(' ')[0]))
muts_g = muts.drop_duplicates(['Patient ID','EFFECT','GENE','POSITION']).copy()
muts_g = muts_g.groupby(['mut_loc','mut_cat']).count().reset_index()

muts_g = pd.DataFrame(index = pd.MultiIndex.from_product([muts_g['mut_loc'].drop_duplicates(),muts_g['mut_cat'].drop_duplicates()] )).merge(muts_g[['mut_loc','mut_cat','Patient ID']], on = ['mut_loc','mut_cat'], how = 'left').fillna(0).reset_index()

muts_g = muts_g.sort_values('mut_loc',ascending = False)

fig,ax = plt.subplots()

bottom = np.zeros(3)

for mc in muts_g['mut_cat'].unique().tolist():
    print(mc)
    tmp = muts_g[muts_g['mut_cat'] == mc].copy()
    ax.bar(tmp['mut_loc'],tmp['Patient ID'], bottom = bottom, color = mut_color_dict.get(mc))
    bottom = bottom + np.array(tmp['Patient ID'].tolist())

ax.set_ylabel('Mutation count')
ax.set_xticklabels(['Shared\nn=44','Primary exclusive\nn=91','Metastatic exclusive\nn=13'])

ax.legend([Patch(color = x) for x in mut_color_dict.values()], list(mut_color_dict.keys()), handlelength = 1)
fig.tight_layout()
fig.savefig("G:/Andy Murtha/Ghent/M1RP/dev/Met_vs_primary/pEx_vs_mEx_counts.pdf")

# =============================================================================
# Boxplot of adjVAF
# =============================================================================

boxplot = [muts[muts['mut_loc'] == ml]['aVAF'].tolist() for ml in muts['mut_loc'].unique().tolist()]
labels = [ml for ml in muts['mut_loc'].unique().tolist()]

fig, ax = plt.subplots()

ax.boxplot(boxplot, labels = labels)

ax.set_ylabel('Clonality (VAF/TC)')
ax.set_xticklabels(['Shared','Primary exclusive','Metastatic exclusive'])

print(stats.kruskal(boxplot[0], boxplot[1], boxplot[2]))
fig.tight_layout()
fig.savefig("G:/Andy Murtha/Ghent/M1RP/dev/Met_vs_primary/pEx_vs_mEx_clonality.pdf")