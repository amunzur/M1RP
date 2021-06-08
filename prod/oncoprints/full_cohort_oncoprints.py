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

genes = ['TP53','PTEN','RB1','FOXA1','SPOP','AR', 'BRCA2','CDK12']
sample_type_color = {0:'green',
                     1:'yellow',
                     2:'purple',
                     3:'red',
                     np.nan:'white'}

tc_color = { 1:'#ffcccc',
             3:'#ff6666',
             5:'#ff0000'}

cn_color = {-2:'#3F60AC',
            -1:'#9CC5E9',
            0:'#E6E7E8',
            1:'#F59496',
            2:'#EE2D24'}

muts_color = {'Missense':'#79B443',
              'Frameshift':'#FFC907',
              'Stopgain':'#FFC907',
              'Non-frameshift':'#BD4398',
              'Splice':'#FFC907'}

def get_tc_color(x):
    if x == 0:
        return 'lightgrey';
    elif x < 0.2:
        return 'grey'
    else:
        return '#ff0000'


n = 2

s1 = ['ID3','ID43','ID10','ID14','ID15','ID19','ID21','ID23','ID33','ID38','ID40', 'ID1', 'ID4', 'ID5','ID6', 'ID7', 'ID9', 'ID16', 'ID17', 'ID18', 'ID20', 'ID22']

s2 = ['ID24', 'ID25', 'ID27', 'ID28', 'ID31', 'ID32', 'ID34', 'ID35', 'ID36', 'ID37', 'ID39', 'ID42', 'ID2', 'ID8', 'ID11', 'ID12', 'ID13', 'ID26', 'ID29', 'ID30', 'ID41']


# =============================================================================
# Import copy number and mutation data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
samples = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

samples.columns = samples.iloc[0]
samples = samples.drop(samples.index[0])

# =============================================================================
# Get all samples. Order by patient then sample type
# =============================================================================

samples['pt_order'] = samples['Patient ID'].map(dict(zip(s1+s2,np.arange(len(s1+s2)))))

samples = samples[~samples['Sample ID'].str.contains('_UCC')]

samples.loc[samples['Sample ID'].str.contains('PB'), 'sample_order'] = '0'
samples.loc[samples['Sample ID'].str.contains('_RP'), 'sample_order'] = '1'
samples.loc[samples['Sample ID'].str.contains('MB'), 'sample_order'] = '2'
samples.loc[samples['Sample ID'].str.contains('MT'), 'sample_order'] = '2'
samples.loc[samples['Sample ID'].str.contains('MLN'), 'sample_order'] = '2'
samples.loc[samples['Sample ID'].str.contains('cfDNA'), 'sample_order'] = '3'
samples['sample_order'] = samples['sample_order'].astype(int)

if n == 1:
    samples = samples[samples['Patient ID'].isin(s1)]
else:
    samples = samples[samples['Patient ID'].isin(s2)]
    

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

# cn['Copy_num'] = cn['Copy_num'].replace({2:1,-2:-1})

# =============================================================================
# Create matrix
# =============================================================================

matrix = pd.DataFrame(index = ['Sample type','Tumor content']+genes, columns = sample_list)


for sample in sample_list:
    if sample == '': continue;
    matrix.at['Sample type', sample] = samples.at[sample, 'sample_order']
    matrix.at['Tumor content', sample] = samples.at[sample,'Final tNGS_TC']
    for gene in genes:
        if gene == 'FOXA1':
            matrix.at[gene,sample] = 0
        else:
            matrix.at[gene,sample] = cn.at[(sample,gene),'Copy_num']
        
matrix = matrix.astype(float)
matrix = matrix.transpose()
matrix = matrix.reset_index()
matrix = matrix.rename(columns = {'index':"Sample ID"})
matrix = matrix[matrix['Sample ID'] != '']

# =============================================================================
# sort matrix
# =============================================================================

matrix['Patient ID'] = matrix['Sample ID'].str.split('_').str[1]
matrix['order'] = matrix['Patient ID'].map(dict(zip(s1+s2,np.arange(len(s1+s2)))))
matrix = matrix.sort_values(['order','Sample type','Tumor content','TP53','RB1','FOXA1','SPOP'], ascending = [True,True,False,True,True,True,True])
matrix['x'] = np.arange(len(matrix))


# =============================================================================
# Separate patients
# =============================================================================

for pt in matrix['order'].unique().tolist():
    matrix.loc[matrix['order'] > pt, 'x'] += 4
    
# =============================================================================
# Get x tick coordinates
# =============================================================================
    
pts = matrix['Patient ID'].unique().tolist()

xticks = matrix[['order','x']].groupby(['order']).mean()['x'].to_list()

matrix = matrix.set_index('x')
matrix = matrix.drop(['Patient ID','order'], axis = 1)

# =============================================================================
# Add x coordinates to mutation data
# =============================================================================

tmp = pd.DataFrame({'x':matrix.index, 'Sample ID':matrix['Sample ID']})
muts = muts.merge(tmp, on = 'Sample ID')
del tmp
muts['EFFECT'] = muts['EFFECT'].str.split(' ').str[0]

# ============================================================================
# Create plot
# =============================================================================

fig,ax = plt.subplots(figsize = (7,1.5))

# =============================================================================
# Plot sample type
# =============================================================================

ax.bar(matrix.index.tolist(), 0.8, width = 1, color = matrix['Sample type'].apply(lambda x: sample_type_color.get(x)).tolist())

# =============================================================================
# Plot gleason
# =============================================================================

ax.bar(matrix.index.tolist(), 0.8, width = 1, bottom = 1, color = matrix['Tumor content'].apply(lambda x: get_tc_color(x)).tolist())

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
ax.set_yticklabels(matrix.columns.tolist()[1:3]+['$\it{%s}$' % s for s in matrix.columns.tolist()[3:]], fontsize = 6)
ax.invert_yaxis()

ax.set_xticks(xticks)
if n == 1:
    ax.set_xticklabels(s1, fontsize = 6, ha = 'center')
elif n == 2:
    ax.set_xticklabels(s2, fontsize = 6, ha = 'center')

# X axis alteration
# ax.set_xlabel('Sample', fontsize = 8)
ax.set_xlim(-0.4, max(matrix.index))
ax.tick_params(labelbottom = True, bottom = False, left = False)
ax.tick_params(axis = 'x', pad = -2)

ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.tight_layout()

fig.subplots_adjust(left = 0.1, right = 0.99)

plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Oncoprints/CIHR_full_oncoprint_%i.pdf' % n)
plt.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/Oncoprints/CIHR_full_oncoprint_%i.png' % n, dpi = 300)