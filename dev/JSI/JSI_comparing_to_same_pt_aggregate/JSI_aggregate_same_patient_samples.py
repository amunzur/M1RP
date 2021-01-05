#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 10:28:26 2020

@author: echen
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns


###################################################################################################################
# Loading input files and cleaning data
###################################################################################################################
baits = pd.read_csv("Data/baits_hg38.bed",sep='\t',names=['CHROM','START','END','GENE'])
genes_on_panel = set(baits['GENE'])
genes_on_panel.add('Intergenic')
mut = pd.read_csv('C:/Users/Emilia Chen/Dropbox/Ghent M1 2019/sandbox/mutations/temporary_clean_data_files/m1rp_mutations.tsv',sep='\t')
# mut = pd.read_csv("Data/m1rp_mutations_new.tsv",sep='\t')

# Restricting to primary samples (i.e. PB and RP)
mut = mut[mut['CLASSIFICATION']=='Somatic']
mut['Sample type'] = mut['Sample ID'].str.split('_').str[2]
mut['Sample type'] = mut['Sample type'].astype(str)
mut = mut[mut['Sample type'].str.contains('PB|RP', regex=True)]
mut['GENE'] = mut['GENE'].astype(str)
mut['GENE'] = np.where(mut['EFFECT']=='Intergenic', 'Intergenic',  mut['GENE'])

# Restricting to genes that are on the panel
for index, row in mut.iterrows():
    for g in mut.at[index,'GENE'].split(';'):
        if g in genes_on_panel:
            mut.at[index,'GENE'] = g
mut = mut[mut['GENE'].isin(genes_on_panel)]

# Restricting to intergenic positions that are on the panel
intergenic = mut[mut['EFFECT']=='Intergenic']
intergenic['onTarget'] = 'F'
for index, row in intergenic.iterrows():
    position = intergenic.at[index,'POSITION']
    chrom = intergenic.at[index,'CHROM']
    tmp = baits[baits['CHROM']==chrom]
    for index_baits, row_baits in tmp.iterrows():
        start = tmp.at[index_baits,'START']
        end = tmp.at[index_baits,'END']
        if start <= position and position <= end:
            intergenic.at[index,'onTarget'] = 'T'
intergenic = intergenic[['Sample ID','EFFECT','CHROM','POSITION','onTarget']]
mut = mut.merge(intergenic,on=['Sample ID','EFFECT','CHROM','POSITION'],how='left')
mut = mut[~(mut['onTarget']=='F')]
mut = mut.reset_index(drop=True)

###################################################################################################################
# Calculating JSI
###################################################################################################################

# Helper
unique_mutation_columns = ['Patient ID','REF','ALT','GENE','EFFECT','CHROM','POSITION']
def get_mutation_jsi_with_pooled(muts, sample):
    patient = sample.split('_')[1]
    muts_tmp = muts[muts['Patient ID']==patient].copy()
    muts_unique = muts_tmp.drop_duplicates(unique_mutation_columns)[unique_mutation_columns]

    muts_unique = muts_unique.set_index(unique_mutation_columns)
    muts_unique['Pooled'] = 0
    muts_unique['Sample'] = 0
    
    muts_tmp = muts_tmp.set_index(unique_mutation_columns)
    muts_sample = muts_tmp[muts_tmp['Sample ID'] == sample]
    
    for index, row in muts_unique.iterrows():
        muts_unique.at[index, 'Pooled'] = 1
        if index in muts_sample.index:
            muts_unique.at[index, 'Sample'] = 1
    muts_unique['Shared'] = 1
    muts_unique.loc[muts_unique['Pooled'] != muts_unique['Sample'], 'Shared'] = 0
    return muts_unique[['Shared']]

all_samples = list(set(mut['Sample ID']))
jsi_sample = pd.DataFrame(columns=['Patient ID','Sample ID','JSI','numMut'])
jsi_sample['Sample ID'] = all_samples
jsi_sample['Patient ID'] = jsi_sample['Sample ID'].str.split('_').str[1]
jsi_sample = jsi_sample.set_index('Sample ID')
for index, row in jsi_sample.iterrows():
    mutation_jsi = get_mutation_jsi_with_pooled(mut, index)
    if len(mutation_jsi) > 0:
        jsi_sample.at[index,'JSI'] = mutation_jsi['Shared'].sum() / len(mutation_jsi)
    else:
        jsi_sample.at[index,'JSI'] = 0
    jsi_sample.at[index,'numMut'] = len(mutation_jsi)


###################################################################################################################
# Ploting results
###################################################################################################################

# Filtering out patients with numMut<=2 and numSample<=2
numMut_threshold = 2
numSample_threshold = 2
jsi_sample_eligible = jsi_sample[jsi_sample['numMut']>numMut_threshold]
numSample = jsi_sample.groupby(by='Patient ID').size().reset_index(name='numSample')
jsi_sample_eligible = jsi_sample_eligible.merge(numSample,on='Patient ID',how='left')
jsi_sample_eligible = jsi_sample_eligible[jsi_sample_eligible['numSample']>numSample_threshold]
jsi_sample_eligible['JSI'] = jsi_sample_eligible['JSI'].astype(float)
median_JSI = jsi_sample_eligible.groupby(by='Patient ID')['JSI'].median().reset_index(name='Median JSI')
jsi_sample_eligible = jsi_sample_eligible.merge(median_JSI,on='Patient ID',how='left')

# Sorting by patient ID for plotting
jsi_sample['indexNumber'] = jsi_sample['Patient ID'].str.replace('ID','')
jsi_sample['indexNumber'] = jsi_sample['indexNumber'].astype(int)
jsi_sample = jsi_sample.sort_values(by=['indexNumber'],ascending=True)

# Ordering
# jsi_sample_eligible = jsi_sample_eligible.sort_values(by='Median JSI',ascending=True)
# jsi_sample_eligible = jsi_sample_eligible.sort_values(by='numSample',ascending=True)


fig = plt.figure(figsize=(8,12))
gs = gridspec.GridSpec(ncols = 2, nrows = 1, width_ratios=[0.4,1],wspace=0.2)

# Ploting JSI swarmplot and boxplot
ax = fig.add_subplot(gs[0,1])
sns.boxplot(y="Patient ID", x="JSI", orient="h", ax=ax,whis=1.5,showfliers = False,boxprops=dict(alpha=.1),linewidth=0.8,zorder=0,color='#999DA0',data=jsi_sample_eligible)
sns.swarmplot(y="Patient ID", x="JSI", orient='h',ax=ax, color='#0a043c',size=3,data=jsi_sample_eligible)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(0,1.1)
ax.set_xlabel('JSI')
ax.set_ylabel('')

# Ploting num unique mut in a patient
ax2 = fig.add_subplot(gs[0,0],sharey=ax)
jsi_sample_eligible['color'] = np.where(jsi_sample_eligible['numMut'] > 2, 'k', '#fecd1a')
ax2.barh(jsi_sample_eligible["Patient ID"],jsi_sample_eligible["numMut"],color=jsi_sample_eligible['color'])
ax2.scatter([],[],color='#fecd1a',label='numMut<=2')
ax2.scatter([],[],color='k',label='numMut>2')
ax2.yaxis.tick_right()
ax2.tick_params(right=True, labelright=False)
ax2.spines['left'].set_visible(False)
ax2.spines['top'].set_visible(False)
plt.gca().invert_xaxis()
ax2.set_xlabel('Num unique muts')
ax2.set_ylabel('')

# fig.savefig('JSI_with_same_pt_aggregate_filtered_sort_by_patientID.png',dpi=600,bbox_inches='tight')