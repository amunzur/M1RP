# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 12:27:16 2021

@author: amurtha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string

mut_genes = ['TP53','FOXA1','BRCA2','SPOP','PTEN','RB1']
cn_genes = ['TP53','PTEN','RB1','BRCA2','NKX3-1','CHD1']

effect_dict = {'Missense':'Missense', 'Stopgain':'Truncating', 'Non-frameshift':'Non-frameshift', 'Frameshift':'Truncating', 'Splice':'Truncating'}

cn_color = {-2:'#3F60AC',
            -1:'#9CC5E9',
            0:'#E6E7E8',
            1:'#F59496',
            2:'#EE2D24'}

mut_colors = {'Missense':'#79B443',
              'Truncating':'#FFC907',
              'Non-frameshift':'#BD4398'}

def keepCodingMutations(df_muts):
    return df_muts[(df_muts['EFFECT'].str.contains("Missense", regex=False)) | (df_muts['EFFECT'].str.contains("Stopgain", regex=False)) | (df_muts['EFFECT'].str.contains("Frameshift", regex=False)) | (df_muts['EFFECT'].str.contains("Splice", regex=False)) | (df_muts['EFFECT'].str.contains("Non-frameshift indel", regex=False)) | (df_muts['EFFECT'] == 'EFFECT') | ((df_muts['EFFECT'] == 'Upstream')&(df_muts['GENE'] == 'TERT'))]

# =============================================================================
# Import data
# =============================================================================

muts = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_mutations_inclDependent.tsv', sep = '\t')
cn = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/copy_number/final melted cna files/M1RP_cna.tsv', sep = '\t')
tc = pd.read_csv('https://docs.google.com/spreadsheets/d/13A4y3NwKhDevY9UF_hA00RWZ_5RMFBVct2RftkSo8lY/export?format=csv&gid=963468022')
gl = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_germline_mutations.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

tc = tc[~tc['Sample ID'].str.contains('cfDNA')]

# =============================================================================
# Keep mutation and CN from TC positive samples
# =============================================================================

tc = tc[tc['Final tNGS_TC'] > 0]
tc_high = tc[tc['Final tNGS_TC'] > 0.2].copy()

muts = muts[muts['Sample ID'].isin(tc['Sample ID'].tolist())]
cn = cn[cn['Sample ID'].isin(tc_high['Sample ID'].tolist())]

# =============================================================================
# Keep mutations and CNA in genes
# =============================================================================

muts = keepCodingMutations(muts)

muts = muts[muts['GENE'].isin(mut_genes)]
cn = cn[cn['GENE'].isin(cn_genes)]
gl = gl[gl['GENE'].isin(mut_genes)]

# =============================================================================
# Get sample counts per patient for CNA
# =============================================================================

high_tc_s_counts = tc_high.groupby('Patient ID').count()[['Sample ID']]
high_tc_s_counts.columns = ['pt_sample_count']

cn = cn.merge(high_tc_s_counts, left_on = 'Patient ID', right_index = True, how = 'left')

# =============================================================================
# Add CNA counts. Keep CNA where count > 3. If not enough samples are present,
# remove it. Including copy neutral
# =============================================================================

cn_counts = cn.groupby(['Patient ID','GENE','Copy_num']).count()[['Log_ratio']]
cn_counts.columns = ['cn_count']
cn_counts = cn_counts.reset_index()
cn_counts = cn_counts[(cn_counts['cn_count'] >= 3)|(cn_counts['Copy_num'].isin([-2,2]))]

cn = cn.merge(cn_counts, on = ['Patient ID','GENE','Copy_num'])

# =============================================================================
# Separate into cn_m and cn_p
# =============================================================================

cn['Sample category'] = cn['Sample ID'].str.split('_').str[2].str.strip(string.digits)
cn['Sample category'] = cn['Sample category'].map({'RP':'Primary', 'MLN':'Metastatic', 'PB':'Primary', 'MT':'Metastatic', 'MB':'Metastatic', 'UCC':'UCC'})

cn = cn[cn['Sample category'] != 'UCC']

cn_m = cn[cn['Sample category'] == 'Metastatic'].copy()
cn_p = cn[cn['Sample category'] == 'Primary'].copy()

# =============================================================================
# Get denominators for both
# =============================================================================

m_total = len(cn_m['Patient ID'].unique())
p_total = len(cn_p['Patient ID'].unique())

# =============================================================================
# Get metastatic cn counts
# =============================================================================

cn_m = cn_m[['Patient ID','GENE','Copy_num']]
cn_m = cn_m[cn_m['Copy_num'] != 0]
cn_m = cn_m.drop_duplicates()

cn_m_dd = cn_m[cn_m['Copy_num'] == -2].copy()
for index, row in cn_m_dd.iterrows():
    cn_m = cn_m[(cn_m['Patient ID'] != row['Patient ID'])|(cn_m['GENE'] != row['GENE'])|(cn_m['Copy_num'] != -1)]
    
cn_m = cn_m.groupby(['GENE','Copy_num']).count().reset_index()
cn_m = cn_m[['GENE','Patient ID','Copy_num']]
cn_m.columns = ['GENE','Count','Copy_num']

cn_m['Total'] = m_total
cn_m['Frequency'] = cn_m['Count'] / cn_m['Total']

# =============================================================================
# Get primary cn counts
# =============================================================================

cn_p = cn_p[['Patient ID','GENE','Copy_num']]
cn_p = cn_p[cn_p['Copy_num'] != 0]
cn_p = cn_p.drop_duplicates()

cn_p_dd = cn_p[cn_p['Copy_num'] == -2].copy()
for index, row in cn_p_dd.iterrows():
    cn_p = cn_p[(cn_p['Patient ID'] != row['Patient ID'])|(cn_p['GENE'] != row['GENE'])|(cn_p['Copy_num'] != -1)]
    
cn_p = cn_p.groupby(['GENE','Copy_num']).count().reset_index()
cn_p = cn_p[['GENE','Patient ID','Copy_num']]
cn_p.columns = ['GENE','Count','Copy_num']

cn_p['Total'] = p_total
cn_p['Frequency'] = cn_p['Count'] / cn_p['Total']

# =============================================================================
# Keep directional
# =============================================================================

cn_p = cn_p[(cn_p['Copy_num'] < 0)|(cn_p['GENE'].isin([]))]
cn_p['color'] = cn_p['Copy_num'].map(cn_color)

cn_m = cn_m[(cn_m['Copy_num'] < 0)|(cn_m['GENE'].isin([]))]
cn_m['color'] = cn_m['Copy_num'].map(cn_color)

# =============================================================================
# Assign x coordinates and start plottin :)
# =============================================================================

x_dict = dict(zip(cn_genes, np.arange(len(mut_genes), len(mut_genes)+len(cn_genes),1)))
cn_p['x'] = cn_p['GENE'].map(x_dict)
cn_m['x'] = cn_m['GENE'].map(x_dict)+0.4

# =============================================================================
# Start plotting
# =============================================================================

fig,ax = plt.subplots(figsize = (3,2.5))

# =============================================================================
# Copy number
# =============================================================================

cn_p_1 = cn_p[cn_p['Copy_num'] == -1]
cn_p_2 = cn_p[cn_p['Copy_num'] == -2]

cn_p_2 = cn_p_2.merge(cn_p_1[['GENE','Frequency']].rename(columns = {'Frequency':'Freq_1'}), on = ['GENE'])

ax.bar(cn_p_1['x'], cn_p_1['Frequency'], color = cn_p_1['color'], width = 0.4)
ax.bar(cn_p_2['x'], cn_p_2['Frequency'], bottom = cn_p_2['Freq_1'], color = cn_p_2['color'], width = 0.4)


cn_m_1 = cn_m[cn_m['Copy_num'] == -1]
cn_m_2 = cn_m[cn_m['Copy_num'] == -2]

cn_m_2 = cn_m_2.merge(cn_m_1[['GENE','Frequency']].rename(columns = {'Frequency':'Freq_1'}), on = ['GENE'])

ax.bar(cn_m_1['x'], cn_m_1['Frequency'], color = cn_m_1['color'], width = 0.4)
ax.bar(cn_m_2['x'], cn_m_2['Frequency'], bottom = cn_m_2['Freq_1'], color = cn_m_2['color'], width = 0.4)

# =============================================================================
# Do mutations
# =============================================================================

# =============================================================================
# Get sample counts per patient for mutations
# =============================================================================

s_counts = tc.groupby('Patient ID').count()[['Sample ID']]
s_counts.columns = ['pt_sample_count']

muts = muts.merge(s_counts, left_on = 'Patient ID', right_index = True, how = 'left')

# =============================================================================
# Add mutation counts. Keep mutations where count > 3 or = num_samples
# =============================================================================

m_counts = muts.groupby(['Patient ID','GENE','POSITION','EFFECT']).count()[['CHROM']]
m_counts.columns = ['Mutation count']
m_counts = m_counts.reset_index()

muts = muts.merge(m_counts, on = ['Patient ID','GENE','POSITION','EFFECT'], how = 'left')

muts = muts[(muts['Mutation count'] >= 3)|((muts['pt_sample_count'] < 3)&(muts['pt_sample_count'] == muts['Mutation count']))]

# =============================================================================
# Split into metastatic and primary
# =============================================================================

tc['Sample category'] = tc['Sample ID'].str.split('_').str[2].str.strip(string.digits)
tc['Sample category'] = tc['Sample category'].map({'RP':'Primary', 'MLN':'Metastatic', 'PB':'Primary', 'MT':'Metastatic', 'MB':'Metastatic', 'UCC':'UCC'})

tc = tc[tc['Sample category'] != 'UCC']

tc_m = tc[tc['Sample category'] == 'Metastatic'].copy()
tc_p = tc[tc['Sample category'] == 'Primary'].copy()

# =============================================================================
# Get denominators for both
# =============================================================================

m_total = len(tc_m['Patient ID'].unique())
p_total = len(tc_p['Patient ID'].unique())

# =============================================================================
# Get counts for metastatic and primary
# =============================================================================

mut_freq_p = muts[muts['Sample ID'].isin(tc_p['Sample ID'])].copy()
mut_freq_p['EFFECT'] = mut_freq_p['EFFECT'].str.split(' ').str[0].map(effect_dict)
mut_freq_p = mut_freq_p[['Patient ID','GENE','EFFECT','CHROM']]
mut_freq_p = mut_freq_p.drop_duplicates(['Patient ID','GENE'])

mut_freq_p = mut_freq_p.groupby(['GENE','EFFECT']).count().reset_index()

mut_freq_p = mut_freq_p[['GENE','EFFECT','Patient ID']]
mut_freq_p.columns = ['GENE','EFFECT','Count']

mut_freq_m = muts[muts['Sample ID'].isin(tc_m['Sample ID'])].copy()
mut_freq_m['EFFECT'] = mut_freq_m['EFFECT'].str.split(' ').str[0].map(effect_dict)
mut_freq_m = mut_freq_m[['Patient ID','GENE','EFFECT','CHROM']]
mut_freq_m = mut_freq_m.drop_duplicates(['Patient ID','GENE'])

mut_freq_m = mut_freq_m.groupby(['GENE','EFFECT']).count().reset_index()

mut_freq_m = mut_freq_m[['GENE','EFFECT','Patient ID']]
mut_freq_m.columns = ['GENE','EFFECT','Count']

# =============================================================================
# Get frequencies
# =============================================================================

mut_freq_p['Total'] = p_total
mut_freq_p['Frequency'] = mut_freq_p['Count'] / mut_freq_p['Total']

mut_freq_m['Total'] = m_total
mut_freq_m['Frequency'] = mut_freq_m['Count'] / mut_freq_m['Total']

# =============================================================================
# Set up color and x coordinates
# =============================================================================

mut_freq_m['color'] = mut_freq_m['EFFECT'].map(mut_colors)
mut_freq_p['color'] = mut_freq_p['EFFECT'].map(mut_colors)


mut_x_dict = dict(zip(mut_genes,np.arange(len(mut_genes))))
mut_freq_p['x'] = mut_freq_p['GENE'].map(mut_x_dict)
mut_freq_m['x'] = mut_freq_m['GENE'].map(mut_x_dict)+0.4

# =============================================================================
# 
# =============================================================================

bottom = pd.DataFrame(index = mut_genes, columns = ['b']).fillna(0)

for effect in list(mut_colors.keys()):
    tmp = mut_freq_p[mut_freq_p['EFFECT'] == effect].copy().set_index(['GENE'])
    tmp = tmp.merge(bottom, left_index = True, right_index = True)
    ax.bar(tmp['x'], tmp['Frequency'], width = 0.4, color = tmp['color'], bottom = tmp['b'])
    bottom = bottom.merge(tmp[['Frequency']], left_index = True, right_index = True, how = 'left').fillna(0)
    bottom['b'] = bottom['b']+bottom['Frequency']
    bottom = bottom.drop('Frequency', axis = 1)
    
bottom = pd.DataFrame(index = mut_genes, columns = ['b']).fillna(0)

for effect in list(mut_colors.keys()):
    tmp = mut_freq_m[mut_freq_m['EFFECT'] == effect].copy().set_index(['GENE'])
    tmp = tmp.merge(bottom, left_index = True, right_index = True)
    ax.bar(tmp['x'], tmp['Frequency'], width = 0.4, color = tmp['color'], bottom = tmp['b'])
    bottom = bottom.merge(tmp[['Frequency']], left_index = True, right_index = True, how = 'left').fillna(0)
    bottom['b'] = bottom['b']+bottom['Frequency']
    bottom = bottom.drop('Frequency', axis = 1)
    
ax.set_xticks(np.arange(0.2,len(mut_genes)+len(cn_genes),1))
ax.set_xticklabels(mut_genes+cn_genes, rotation = 90, fontsize = 6)

ax.set_ylabel('Alteration frequency', fontsize = 6)
ax.set_ylim(0,0.7)

ax.tick_params(labelsize = 6, bottom = False)

fig.tight_layout()
fig.subplots_adjust(left = 0.13, right = 0.98, top = 0.96, bottom = 0.22)

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/summary/AlterationFrequencyByMetVsPri.pdf')