# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 11:15:04 2020
@author: amurtha
"""

import pandas as pd
import itertools
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy

# =============================================================================
# Constants
# =============================================================================

cohort = 'M1RP'
unique_mutation_columns = ['Patient ID','REF','ALT','GENE','EFFECT','CHROM','POSITION']
patient = 'ID1'

# =============================================================================
# Helpers
# =============================================================================

def get_mutation_jsi(muts, samples):
    muts_tmp = muts[muts['Sample ID'].isin(samples)].copy()
    muts_unique = muts_tmp.drop_duplicates(unique_mutation_columns)[unique_mutation_columns]
    muts_unique = muts_unique.set_index(unique_mutation_columns)
    muts_unique['Sample1'] = 0
    muts_unique['Sample2'] = 0
    muts_tmp = muts_tmp.set_index(unique_mutation_columns)
    muts_sample1 = muts_tmp[muts_tmp['Sample ID'] == samples[0]]
    muts_sample2 = muts_tmp[muts_tmp['Sample ID'] == samples[1]]
    for index, row in muts_unique.iterrows():
        if index in muts_sample1.index:
            muts_unique.at[index, 'Sample1'] = 1
        if index in muts_sample2.index:
            muts_unique.at[index, 'Sample2'] = 1
    muts_unique['Shared'] = 1
    muts_unique.loc[muts_unique['Sample1'] != muts_unique['Sample2'], 'Shared'] = 0
    return muts_unique[['Shared']]

def get_cn_jsi(cn,samples):
    cn_tmp = cn[cn['Sample ID'].isin(samples)].copy()
    cn_tmp['Copy_num'] = cn_tmp['Copy_num'].replace({2:1,-2:-1})
    cn_tmp = cn_tmp.pivot(index = 'GENE', columns = 'Sample ID', values = 'Copy_num')
    cn_tmp = cn_tmp[(cn_tmp[samples[0]] != 0)|(cn_tmp[samples[1]] != 0)]
    cn_tmp['Shared'] = 1
    cn_tmp.loc[cn_tmp[samples[0]] != cn_tmp[samples[1]], 'Shared'] = 0
    return cn_tmp[['Shared']]

# =============================================================================
# Import mutations, CN, TC, and clinical data
# =============================================================================

muts = pd.read_csv('C:/Users/Emilia Chen/Dropbox/Ghent M1 2019/sandbox/mutations/final melted mutations/%s_mutations_inclDependent.tsv' % cohort, sep = '\t')
cn = pd.read_csv('C:/Users/Emilia Chen/Dropbox/Ghent M1 2019/sandbox/copy number/final melted cna files/%s_allSamples_cna.tsv' % cohort, sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc = tc[tc['Cohort'] == cohort]
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].str.split('%').str[0].astype(float) / 100
tc = tc[tc['Final tNGS_TC'] > 0]

# Test run on ID1
tc = tc[tc['Patient ID'].isin([patient])]

tc_dict = dict(zip(tc['Sample ID'].tolist(), tc['Final tNGS_TC'].tolist()))
# =============================================================================
# Keep samples with TC > 40%
# =============================================================================

# tc = tc[tc['Final tNGS_TC'] > 0.4]

# =============================================================================
# Create same patient sample combination list
# =============================================================================

ls = []

# test run on ID1
for pt in [patient]:
    samples = tc[tc['Patient ID'] == pt]['Sample ID'].tolist()
    ls.append(list(itertools.product(samples, samples)))
    
samples = []
for l in ls:
    samples = samples+l
    
del pt, ls

samples = [sample for sample in samples if sample[0] != sample[1]]
index = pd.MultiIndex.from_tuples(samples, names=['Sample1', 'Sample2'])

# =============================================================================
# Create dataframe
# =============================================================================

jsi = pd.DataFrame(index = index)
matrix = pd.DataFrame(columns = tc['Sample ID'].unique(),index = tc['Sample ID'].unique())

for index, row in jsi.iterrows():
    sample1 = index[0]
    sample2 = index[1]
    mutation_jsi = get_mutation_jsi(muts, [sample1,sample2])
    cn_jsi = get_cn_jsi(cn,[sample1,sample2])
    if len(cn_jsi)+ len(mutation_jsi) > 0:
        if tc_dict.get(sample1) > 0.2 and tc_dict.get(sample2) > 0.2:
            value = (cn_jsi['Shared'].sum()+mutation_jsi['Shared'].sum())/(len(cn_jsi)+ len(mutation_jsi))
        else:
            value = mutation_jsi['Shared'].sum() / len(mutation_jsi)
    else:
        value = 0
    matrix.at[sample1, sample2] = value
# =============================================================================
# Create and plot per patient matrix
# =============================================================================
#  for pt in tc['Patient ID'].unique()
# Test run on ID1
for pt in [patient]:
    samples = tc[tc['Patient ID'] == pt]['Sample ID'].tolist()
    fig,[[ax1,ax2],[ax3,ax4]] = plt.subplots(ncols = 2, nrows = 2, sharey = 'row', gridspec_kw = {'width_ratios':[4,1], 'height_ratios':[1,4],'wspace':0.01, 'hspace':0.1}, figsize = (10,10))
    

#     matrix_tmp = matrix[matrix.index.isin(samples)][samples].fillna(1).to_numpy()
#     matrix_tmp = linkage(matrix_tmp)
#     dend = dendrogram(matrix_tmp, ax = ax1)
#     ax1.set_xticks(np.arange(0,len(samples),1))
#     order = order.sort_values('order',ascending = True).reset_index().drop('index',axis=1)

    # ploting dendrogram
    matrix_tmp = matrix[matrix.index.isin(samples)][samples].fillna(1)
    Z = hierarchy.linkage(matrix_tmp, 'ward')
    dend = hierarchy.dendrogram(Z, ax = ax1, leaf_rotation=90, leaf_font_size=8, labels=matrix_tmp.index)


    # ploting JSI matrix
    order = pd.DataFrame({'order':dend.get('ivl')})

    # ax1.tick_params(bottom = False, labelbottom = False)
    samples = order['order'].tolist()
    
    for x in np.arange(0,len(samples),1):
        for y in np.arange(x,len(samples),1):
            sample1 = samples[x]
            sample2 = samples[y]
            if sample1 == sample2:
                val = 0
            else:
                val = 1 - matrix.at[sample1, sample2]
            ax3.bar(x,0.8, bottom = y, color = str(round(val, 3)))
    ax3.set_xticks(np.arange(0,len(samples),1))
    ax3.set_xticklabels(samples, rotation = 90)
    ax3.set_yticks(np.arange(0.4,len(samples),1))
    ax3.set_yticklabels(samples)
    
    # ploting TC
    tc = tc.set_index('Sample ID')
    tc = tc.reindex(samples)
    tc = tc.reset_index()
    
    ypos =np.arange(0.4,len(samples),1)
    ax4.barh(ypos, tc['Final tNGS_TC'],color = '#808080')
    ax4.set_xlim(0,1)
    ax4.set_xlabel('TC%', fontsize=12)
    ax4.set_xticks([0, 0.20, 0.40, 0.60, 0.80, 1])
    ax4.set_xticklabels([0, 20, 40, 60, 80, 100],fontsize=12)

    fig.savefig(patient+'_jsi.pdf',bbox_inches='tight')