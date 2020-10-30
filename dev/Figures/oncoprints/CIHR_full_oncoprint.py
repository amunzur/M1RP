# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 12:20:35 2020

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# helpers
# =============================================================================

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]

# =============================================================================
# Constants
# =============================================================================

genes = ['TP53','PTEN','RB1','FOXA1','SPOP','AR']
sample_type_color = {0:'green',
                     1:'yellow',
                     2:'purple',
                     3:'red',
                     np.nan:'white'}

ggg_color = {0:'#ffffff',
             1:'#ffcccc',
             2:'#ff9999',
             3:'#ff6666',
             4:'#ff3232',
             5:'#ff0000'}

cn_color = {-2:'#3F60AC',
            -1:'#9CC5E9',
            0:'#E6E7E8',
            1:'#F59496',
            2:'#EE2D24'}

muts_color = {'Missense':'#79B443',
              'Frameshift':'#FFC907',
              'Stopgain':'#FFC907',
              'Non-frameshift':'#BD4398'}

# =============================================================================
# Import copy number and mutation data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/M1RP_allSamples_cna.tsv', sep = '\t')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=0')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

# =============================================================================
# Get all samples. Order by patient then sample type
# =============================================================================

samples = samples[samples['Cohort'] == 'M1RP']
samples = samples[samples['Targeted sequencing'] == 'Completed']
samples = samples[samples['Sample type'].isin(['FFPE tissue', 'cfDNA'])]
samples['pt_order'] = samples['Patient ID'].str.split('ID').str[1].astype(int)

samples.loc[samples['Sample ID'].str.contains('PB'), 'sample_order'] = '0'
samples.loc[samples['Sample ID'].str.contains('_RP'), 'sample_order'] = '1'
samples.loc[samples['Sample ID'].str.contains('MB'), 'sample_order'] = '2'
samples.loc[samples['Sample ID'].str.contains('MLN'), 'sample_order'] = '2'
samples.loc[samples['Sample ID'].str.contains('cfDNA'), 'sample_order'] = '3'
samples['sample_order'] = samples['sample_order'].astype(int)

samples = samples[samples['pt_order'] < 24]
samples = samples[samples['pt_order'] != 20]
samples = samples.sort_values(['pt_order','sample_order']).reset_index(drop = True)

# =============================================================================
# Add space between samples of different patients
# =============================================================================

sample_list = []
for index,sample in samples['Sample ID'].iteritems():
    pt = sample.split('_')[1]
    sample_list.append(sample)
    if index < len(samples)-1 and pt != samples.at[index+1, 'Patient ID']:
        sample_list = sample_list + [''] * 5
        
samples = samples.set_index('Sample ID')

# =============================================================================
# Keep mutation and copy number alterations in genes
# =============================================================================

cn = cn[cn['GENE'].isin(genes)]
muts = muts[muts['GENE'].isin(genes)]
muts = keepCodingMutations(muts)

cn = cn.set_index(['Sample ID','GENE'])

# =============================================================================
# Create matrix
# =============================================================================

matrix = pd.DataFrame(index = ['Sample type','GGG']+genes, columns = sample_list)


for sample in sample_list:
    if sample == '': continue;
    matrix.at['Sample type', sample] = samples.at[sample, 'sample_order']
    matrix.at['GGG',sample] = samples.at[sample,'GGG']
    for gene in genes:
        matrix.at[gene,sample] = cn.at[(sample,gene),'Copy_num']
        
matrix = matrix.astype(float)
matrix = matrix.transpose()
matrix = matrix.reset_index()
matrix = matrix.rename(columns = {'index':"Sample ID"})
matrix = matrix[matrix['Sample ID'] != '']
matrix['GGG'] = matrix['GGG'].fillna(0)
matrix.loc[matrix['Sample type'] > 1, 'GGG'] = 0

# =============================================================================
# Organize mutation
# =============================================================================

tmp = pd.DataFrame({'x':matrix.index, 'Sample ID':matrix['Sample ID']})
muts = muts.merge(tmp, on = 'Sample ID')
del tmp
muts['EFFECT'] = muts['EFFECT'].str.split(' ').str[0]

# =============================================================================
# sort matrix
# =============================================================================

matrix = matrix.sort_values(['Sample type','GGG','TP53','RB1','FOXA1','SPOP'], ascending = [True, False, True,True,True,True])

# ============================================================================
# Create plot
# =============================================================================

fig,ax = plt.subplots(figsize = (7.5,1.25))

# =============================================================================
# Plot sample type
# =============================================================================

ax.bar(matrix.index.tolist(), 0.8, width = 1, color = matrix['Sample type'].apply(lambda x: sample_type_color.get(x)).tolist())

# =============================================================================
# Plot gleason
# =============================================================================

ax.bar(matrix.index.tolist(), 0.8, width = 1, bottom = 1, color = matrix['GGG'].apply(lambda x: ggg_color.get(x)).tolist())

# =============================================================================
# Plot copy number
# =============================================================================

for y, gene in enumerate(genes):
    ax.bar(matrix.index.tolist(), 0.8, width = 1, bottom = 2 + y, color = matrix[gene].apply(lambda x: cn_color.get(x)).tolist())
    
    # =============================================================================
    # Plot mutations
    # =============================================================================
    tmp = muts[muts['GENE'] == gene]
    ax.bar(tmp['x'].tolist(), 0.2, width = 1, bottom = 2.3 + y, color = tmp['EFFECT'].apply(lambda x: muts_color.get(x)).tolist())    

# =============================================================================
# Aethetics
# =============================================================================
    
# Y axis alteration
ax.set_yticks(np.arange(.4, len(matrix.columns)-1, 1))
ax.set_ylim(0,len(matrix.columns)-1)
ax.set_yticklabels(matrix.columns.tolist()[1:3]+['$\it{%s}$' % s for s in matrix.columns.tolist()[3:]], fontsize = 8)
ax.invert_yaxis()

# X axis alteration
# ax.set_xlabel('Sample', fontsize = 8)
ax.set_xlim(-0.4, max(matrix.index))
ax.tick_params(labelbottom = False, bottom = False, left = False)

ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.tight_layout()

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/CIHR_full_oncoprint.pdf')