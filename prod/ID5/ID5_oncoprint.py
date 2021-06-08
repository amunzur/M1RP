# -*- coding: utf-8 -*-
"""
Created on Fri May  7 14:03:37 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import string
import numpy as np

# =============================================================================
# Helpers
# =============================================================================

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False))]

# =============================================================================
# Constants
# =============================================================================

sample_cat = {'RP':1,
              'PB':1,
              'MLN': 2,
              'MB':2,
              'cfDNA': 3}

mut_effect_dict = {'Missense':'#79B443',
              'Frameshift':'#FFC907',
              'Stopgain':'#FFC907',
              'Non-frameshift':'#BD4398',
              'Splice':'#FFC907'}
'''
order = ['M1RP_ID5_RP1',
 'M1RP_ID5_RP2',
 'M1RP_ID5_RP3',
 'M1RP_ID5_RP4',
 'M1RP_ID5_RP5',
 'M1RP_ID5_RP6',
 'M1RP_ID5_RP7',
 'M1RP_ID5_RP8',
 'M1RP_ID5_RP9',
 'M1RP_ID5_RP10',
 'M1RP_ID5_PB1',
 'M1RP_ID5_PB2',
 'M1RP_ID5_PB4',
 'M1RP_ID5_MLN1',
 'M1RP_ID5_MLN2',
 'M1RP_ID5_MLN4',
 'M1RP_ID5_MLN5',
 'M1RP_ID5_MLN6',
 'M1RP_ID5_MLN7',
 'M1RP_ID5_cfDNA_2018Mar30']'''

cn_color_dict = {-2:'#3F60AC',-1:'#9CC5E9', 0:'#E6E7E8', 1:'#F59496', 2:'#EE2D24'}

# =============================================================================
# Import data
# =============================================================================

tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
gl = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_germline_mutations.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)
tc = tc[tc['Final tNGS_TC'] >= 0.20]

tc['Sample Category'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits).apply(lambda x: sample_cat.get(x))

# =============================================================================
# Limit to ID5 TC > 0.2
# =============================================================================

tc = tc[tc['Final tNGS_TC'] > 0.2]
tc = tc[tc['Patient ID'] == 'ID5']
muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]
cn = cn[cn['Sample ID'].isin(tc['Sample ID'].tolist())]

# =============================================================================
# Combine germline and somatic mutations
# =============================================================================

gl = gl[gl['Patient ID'] == 'ID5']
gl = gl[gl['NOTES'].str.contains('ClinVar:Pathogenic')]

muts['Somatic'] = True

for index, row in gl.iterrows():
    muts = muts.append(pd.DataFrame({'Cohort':'M1RP', 
                   'Patient ID':'ID5', 
                   'Sample ID':tc['Sample ID'],
                   'CHROM':row['CHROM'], 
                   'POSITION':row['POSITION'], 
                   'REF':row['REF'], 
                   'ALT':row['ALT'],
                   'GENE':row['GENE'], 
                   'EFFECT':row['EFFECT'], 
                   'Allele_frequency':row['Allele_frequency'], 
                   'Read_depth':row['Read_depth'], 
                   'NOTES':row['NOTES'],
                   'Independent':True, 
                   'Final tNGS_TC':tc['Final tNGS_TC'], 
                   'Clonal':True,
                   'Somatic':False}), ignore_index = True)

# =============================================================================
# Compute the X coordinate for the mutaiton
# =============================================================================

# tc['order'] = tc['Sample ID'].map(dict(zip(order,np.arange(len(order)))))

tc = tc.sort_values(['Sample Category','Final tNGS_TC'], ascending = [True,False])
s_order = tc['Sample ID'].tolist()

muts['x'] = muts['Sample ID'].map(dict(zip(s_order,np.arange(0.,len(tc),1))))
cn['x'] = cn['Sample ID'].map(dict(zip(s_order,np.arange(0.,len(tc),1))))
tc['x'] = tc['Sample ID'].map(dict(zip(s_order,np.arange(0.,len(tc),1))))

# =============================================================================
# Compute Y coordinates
# =============================================================================

genes = ['TP53','PTEN','RB1','KMT2C','SMARCA1']

muts = keepCodingMutations(muts)
muts = muts[muts['GENE'].isin(genes)]
cn = cn[cn['GENE'].isin(genes)]

muts['y'] = muts['GENE'].map(dict(zip(genes, np.arange(0.4,len(genes),1))))
cn['y'] = cn['GENE'].map(dict(zip(genes, np.arange(0.,len(genes),1.))))

# =============================================================================
# Find double mutations
# =============================================================================

mut_count = muts[['GENE','EFFECT']].drop_duplicates().copy()
mut_count = mut_count.groupby(['GENE']).count()
mut_count.columns = ['count']

offset = 0.15

for index, row in mut_count.iterrows():
    if row['count'] > 1:
        effects = muts[muts['GENE'] == index]['EFFECT'].drop_duplicates().to_list()
        for i, effect in enumerate(effects):
            if i == 0:
                muts.loc[(muts['GENE'] == index) & (muts['EFFECT'] == effect), 'x'] = muts['x']+offset
                muts.loc[(muts['GENE'] == index) & (muts['EFFECT'] == effect), 'y'] = muts['y']+offset
            else:
                muts.loc[(muts['GENE'] == index) & (muts['EFFECT'] == effect), 'x'] = muts['x']-offset
                muts.loc[(muts['GENE'] == index) & (muts['EFFECT'] == effect), 'y'] = muts['y']-offset
                
# =============================================================================
# Label colors
# =============================================================================
        
muts['color'] = muts['EFFECT'].str.split(' ').str[0].map(mut_effect_dict)
cn['color'] = cn['Copy_num'].map(cn_color_dict)

# =============================================================================
# Create plot
# =============================================================================
        
fig,[ax2,ax1] = plt.subplots(nrows = 2, figsize = (3,2.25), sharex = True, gridspec_kw={'height_ratios':[0.25,1]})

# =============================================================================
# Plot mutations
# =============================================================================

for b in [True, False]:
    tmp = muts[muts['Somatic'] == b]
    if b == True:
        ax1.scatter(tmp['x'],tmp['y'],c=tmp['color'],marker='s',lw=0, zorder = 100, s = 20)
    else:
        ax1.scatter(tmp['x'],tmp['y'],c=tmp['color'],marker='*',lw=0, zorder = 100, s = 30) 

# =============================================================================
# Plot CNA
# =============================================================================

ax1.bar(cn['x'], 0.8, bottom = cn['y'], color = cn['color'], zorder = 10)

# =============================================================================
# axes adjustments
# =============================================================================

ax1.set_yticks(np.arange(0.4,len(genes),1))
ax1.set_yticklabels(genes, fontsize = 6)

ax1.set_xticks(np.arange(0,len(tc),1))
ax1.set_xticklabels(pd.Series(s_order).str.split('_').str[2], rotation = 90, fontsize = 6)

ax1.spines['left'].set_visible(False)
ax1.spines['bottom'].set_visible(False)

ax1.set_xlim(-0.5,len(s_order)-0.5)
ax1.set_ylim(-0.1,len(genes))

ax1.tick_params(left = False, bottom = False, pad = 0)

ax1.invert_yaxis()

# =============================================================================
# Plot tumor content
# =============================================================================

ax2.bar(tc['x'],tc['Final tNGS_TC'], zorder = 100)

ax2.set_ylim(0,1)
ax2.set_yticks(np.arange(0,1.01,0.25))
ax2.tick_params(labelsize = 6)
ax2.grid(axis = 'y', lw = 0.5, ls = 'dashed', zorder = 0)

ax2.set_ylabel('Tumor content', fontsize = 6)
ax2.spines['bottom'].set_zorder(1000)

fig.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Patient specific/ID5/oncoprint.pdf')