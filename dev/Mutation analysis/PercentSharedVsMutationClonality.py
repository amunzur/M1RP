# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 19:02:03 2021

@author: amurtha
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 14:52:38 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import string
import seaborn as sns

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
                'Non-frameshift':'Non-framshift',
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

# tc['Sample type'] = tc['Sample ID'].apply(lambda x: s_type_dict.get(x.split('_')[2].rstrip(string.digits)))
# tc = tc[tc['Sample type'] == 'Metastatic']


muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]

# =============================================================================
# Create full mutation dataframe. This means that for every TC positive sample
# there is a row for each mutation found in patient pt's samples 
# =============================================================================

# For syntax reasons, intergeneic mutaitons with nan gene must be changed to a value
muts['GENE'] = muts['GENE'].fillna('intergenic')
muts = muts.merge(tc[['Sample ID','Final tNGS_TC']], on = 'Sample ID', how = 'left')

# =============================================================================
# Group by patient, sample type, and mutation. Count frequency.
# =============================================================================

mut_counts = muts.groupby(['Patient ID','GENE','EFFECT','POSITION']).count().reset_index()
mut_counts = mut_counts[['Patient ID','GENE','EFFECT','POSITION','Sample ID']].rename(columns = {'Sample ID':'mut_count'})

mut_counts['Total_samples'] = mut_counts['Patient ID'].apply(lambda x: len(tc[(tc['Patient ID'] == x)&(tc['Cohort'] == 'M1RP')]))

mut_counts['Percent_mut'] = mut_counts['mut_count'] / mut_counts['Total_samples']

mut_counts = mut_counts[mut_counts['Total_samples'] > 5]

# =============================================================================
# Group by patient and mutation. Average clonality
# =============================================================================

muts['aVAF'] = muts['Allele_frequency'] / muts['Final tNGS_TC']

mut_clonality_std =  muts.groupby(['Patient ID','GENE','EFFECT','POSITION']).std().reset_index()
mut_clonality_mean =  muts.groupby(['Patient ID','GENE','EFFECT','POSITION']).mean().reset_index()

mut_clonality_mean['aVAF.mean'] = mut_clonality_mean['aVAF']
mut_clonality_std['aVAF.std'] = mut_clonality_std['aVAF']

mut_counts = mut_counts.merge(mut_clonality_std[['Patient ID','GENE','EFFECT','POSITION','aVAF.std']], on =['Patient ID','GENE', 'EFFECT', 'POSITION'], how = 'left')
mut_counts = mut_counts.merge(mut_clonality_mean[['Patient ID','GENE','EFFECT','POSITION','aVAF.mean']], on =['Patient ID','GENE', 'EFFECT', 'POSITION'], how = 'left')

fig,[ax1,ax2] = plt.subplots(nrows = 2, figsize = (4,5), sharex = True)

ax1.scatter(mut_counts['Percent_mut'], mut_counts['aVAF.mean'], s = 8, alpha = 0.8, color = 'k')
ax2.scatter(mut_counts['Percent_mut'], mut_counts['aVAF.std'], s = 8, alpha = 0.8, color = 'k')

ax1.plot([-1,2],[0.25,0.25], ls = '--', lw = 0.8)

ax1.set_xlim(-0.03, 1.02)

ax1.set_ylabel('Average clonality')

ax2.set_ylabel('Clonality deviation')
ax2.set_xlabel('Fraction of TC+ samples with mutation')

mut_counts = mut_counts[~mut_counts['aVAF.std'].isnull()]
s = stats.linregress(mut_counts['Percent_mut'],mut_counts['aVAF.std'])

fig.align_ylabels()
fig.tight_layout()
fig.savefig("G:/Andy Murtha/Ghent/M1RP/dev/Met_vs_primary/PercentShared_vs_Clonality.pdf")