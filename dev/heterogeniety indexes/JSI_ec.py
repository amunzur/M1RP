# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 11:15:04 2020
@author: amurtha
modified by echen
"""

import pandas as pd
import itertools
import matplotlib as mpl
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

samples = tc[tc['Patient ID'] == patient]['Sample ID'].tolist()
matrix_tmp = matrix[matrix.index.isin(samples)][samples].fillna(1)

# =============================================================================
# Removing M1RP_IDX from columns and rows for clarity
# =============================================================================
for i, col in enumerate(matrix_tmp.columns):
    if 'cfDNA' not in col:
        col = col.split("_")[2]
    else:
        col = col.split("_")[2] + col.split("_")[3]
    matrix_tmp = matrix_tmp.rename(columns={matrix_tmp.columns[i]: col})
    
matrix_tmp = matrix_tmp.reset_index()
for index, row in matrix_tmp.iterrows():
    name = matrix_tmp.at[index,'index']
    if 'cfDNA' not in name:
        matrix_tmp.at[index,'index'] = name.split("_")[2]
    else:
        matrix_tmp.at[index,'index'] = name.split("_")[2] + name.split("_")[3]
matrix_tmp = matrix_tmp.set_index('index')

for index, row in tc.iterrows():
    name = tc.at[index,'Sample ID']
    if 'cfDNA' not in name:
        tc.at[index,'Sample ID'] = name.split("_")[2]
    else:
        tc.at[index,'Sample ID'] = name.split("_")[2] + name.split("_")[3]
# =============================================================================
# Create and plot per patient matrix
# =============================================================================

g = sns.clustermap(matrix_tmp,
                   figsize=(18, 9),
                   metric="correlation",
                   cmap="BuPu",
                   cbar_pos=None)
g.ax_row_dendrogram.set_visible(False)
g.gs.update(left=0.03, right=0.48)
ax = g.ax_heatmap
ax.set_ylabel("")
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)


gs2 = matplotlib.gridspec.GridSpec(1,3, left=0.54,top=0.795) 
ax2 = g.fig.add_subplot(gs2[0])
ax2.set_xlim((0,1))
ax2.tick_params(left = False, labelleft = False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

samples = list()
a = g.ax_heatmap.get_yticklabels()
for x in a:
    samples.append(x.get_text())
tc = tc.set_index('Sample ID')
tc = tc.reindex(samples)
tc = tc.reset_index()

sns.barplot(x="Final tNGS_TC",y='Sample ID',data=tc, orient="h", color='grey', ax = ax2,ci=None)
ax2.set(xlabel="TC", ylabel = "")
g.fig.suptitle(cohort+"_"+patient, x=0.4,y=1,fontsize=16)



filename = cohort + "_" + patient + "_JSI.pdf"
g.savefig(filename)