# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 11:36:35 2021

@author: amurtha
"""

import pandas as pd
import itertools as it
import numpy as np
import matplotlib.pyplot as plt

mut_genes = ['TP53','FOXA1','BRCA2','SPOP','PTEN','RB1']
cn_genes = ['TP53','PTEN','RB1','BRCA2','NKX3-1']

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
clin = pd.read_excel("C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/clinical/ClinicalDataWithLatitude.xlsx")
gl = pd.read_csv('C:/Users/amurtha/Dropbox/Ghent M1 2019/Mar2021_datafreeze/mutations/final melted mutations/M1RP_germline_mutations.tsv', sep = '\t')

tc.columns = tc.iloc[0]
tc = tc.drop(tc.index[0])
tc['Final tNGS_TC'] = tc['Final tNGS_TC'].astype(float)

tc = tc[~tc['Sample ID'].str.contains('cfDNA')]

# =============================================================================
# Clinical
# =============================================================================

clin['study_id'] = clin['study_id'].str.split('-').str[1]
clin = clin[clin['cohort'] == 'M1RP']
clin = clin[['study_id','latitude']]
clin['color'] = clin['latitude'].map({'Low-risk':'#E6E7E8','High-risk':'k'})

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
# Add germline mutations
# =============================================================================
'''
gl = gl[(gl['NOTES'].str.contains('ClinVar:Pathogenic')) | (gl['EFFECT'].isin(['Frameshift','Stopgain','Splice']))]
gl = gl[~gl['NOTES'].str.contains('Benign')]

muts['Somatic'] = True

for index, row in gl.iterrows():
    muts = muts.append(pd.DataFrame({'Cohort':'M1RP', 
                   'Patient ID':row['Patient ID'], 
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
'''
# =============================================================================
# Get sample counts per patient for mutations
# =============================================================================

s_counts = tc.groupby('Patient ID').count()[['Sample ID']]
s_counts.columns = ['pt_sample_count']

muts = muts.merge(s_counts, left_on = 'Patient ID', right_index = True, how = 'left')

# =============================================================================
# Get sample counts per patient for CNA
# =============================================================================

high_tc_s_counts = tc_high.groupby('Patient ID').count()[['Sample ID']]
high_tc_s_counts.columns = ['pt_sample_count']

cn = cn.merge(high_tc_s_counts, left_on = 'Patient ID', right_index = True, how = 'left')

# =============================================================================
# Add mutation counts. Keep mutations where count > 3 or = num_samples
# =============================================================================

m_counts = muts.groupby(['Patient ID','GENE','POSITION','EFFECT']).count()[['CHROM']]
m_counts.columns = ['Mutation count']
m_counts = m_counts.reset_index()

muts = muts.merge(m_counts, on = ['Patient ID','GENE','POSITION','EFFECT'], how = 'left')

muts = muts[(muts['Mutation count'] >= 3)|((muts['pt_sample_count'] < 3)&(muts['pt_sample_count'] == muts['Mutation count']))]

# =============================================================================
# Add CNA counts. Keep CNA where count > 3. If not enough samples are present,
# remove it. Including copy neutral
# =============================================================================

cn_counts = cn.groupby(['Patient ID','GENE','Copy_num']).count()[['Log_ratio']]
cn_counts.columns = ['cn_count']
cn_counts = cn_counts.reset_index()
cn_counts = cn_counts[(cn_counts['cn_count'] >= 3)|(cn_counts['Copy_num'].isin([-2,2]))]

cn = cn.merge(cn_counts, on = ['Patient ID','GENE','Copy_num'])

cn_counts = cn_counts.merge(clin[['study_id','latitude']], left_on = 'Patient ID', right_on = 'study_id')
cn_risk_totals = cn_counts.copy()
cn_risk_totals = cn_risk_totals[['Patient ID','latitude']].drop_duplicates()
cn_risk_totals = cn_risk_totals.groupby('latitude').count()
cn_risk_totals = dict(zip(cn_risk_totals.index, cn_risk_totals['Patient ID']))

## If using 3 as limit for AMP/DeepDel, then use this code so no shallow dels are lost
'''
for index, row in cn.iterrows():
    pt = row['Patient ID']
    g = row['GENE']
    if row['Copy_num'] == -2 and row['cn_count'] < 3:
        cn.loc[(cn['Patient ID'] == pt)&(cn['GENE'] == g)&(cn['Copy_num'] == -1), 'cn_count'] = cn['cn_count']+row['cn_count']
    elif row['Copy_num'] == 2 and row['cn_count'] < 3:
        cn.loc[(cn['Patient ID'] == pt)&(cn['GENE'] == g)&(cn['Copy_num'] == 1), 'cn_count'] = cn['cn_count'] + row['cn_count']
cn = cn[cn['cn_count'] >= 3]
'''
# Keep amps and deep dels and gains/shallow losses if count >= 3
cn = cn[(cn['Copy_num'].isin([-2,2]))|(cn['cn_count'] >= 3)]

# =============================================================================
# Get mut frequencies
# =============================================================================

mut_freq = muts.copy()
mut_freq['EFFECT'] = mut_freq['EFFECT'].str.split(' ').str[0].map(effect_dict)
mut_freq = mut_freq[['Patient ID','GENE','EFFECT','CHROM']]
mut_freq = mut_freq.drop_duplicates(['Patient ID','GENE'])

mut_freq = mut_freq.merge(clin[['study_id','latitude']], left_on = 'Patient ID',right_on = 'study_id', how = 'left').drop('study_id', axis = 1)

mut_freq = mut_freq.groupby(['GENE','EFFECT','latitude']).count().reset_index()

mut_freq = mut_freq[['GENE','EFFECT','latitude','Patient ID']]
mut_freq.columns = ['GENE','EFFECT','latitude','Count']

# =============================================================================
# Get cna freq
# =============================================================================

cn_freq = cn.copy()
cn_freq = cn_freq[['Patient ID','GENE','Copy_num']]
cn_freq = cn_freq[cn_freq['Copy_num'] != 0]
cn_freq = cn_freq.drop_duplicates()

cn_freq_dd = cn_freq[cn_freq['Copy_num'] == -2].copy()
for index, row in cn_freq_dd.iterrows():
    cn_freq = cn_freq[(cn_freq['Patient ID'] != row['Patient ID'])|(cn_freq['GENE'] != row['GENE'])|(cn_freq['Copy_num'] != -1)]

cn_freq = cn_freq.merge(clin[['study_id','latitude']], left_on = 'Patient ID', right_on = 'study_id', how = 'left').drop('study_id', axis = 1)
cn_freq = cn_freq.groupby(['GENE','latitude','Copy_num']).count().reset_index()


cn_freq = cn_freq[['GENE','latitude','Patient ID','Copy_num']]
cn_freq.columns = ['GENE','latitude','Count','Copy_num']

# =============================================================================
# Set x-coordinates
# =============================================================================

mut_x_dict = dict(zip(mut_genes,np.arange(len(mut_genes))))

mut_freq['x'] = mut_freq['GENE'].map(mut_x_dict)
mut_freq.loc[mut_freq['latitude'] == 'High-risk', 'x'] = mut_freq['x']+0.4

# =============================================================================
# Adjust to frequency rather than count
# =============================================================================

mut_risk_totals = {'High-risk':len(clin[clin['latitude'] == 'High-risk']),'Low-risk':len(clin[clin['latitude'] == 'Low-risk'])}
mut_freq['Total'] = mut_freq['latitude'].map(mut_risk_totals)
mut_freq['Frequency'] = mut_freq['Count']/mut_freq['Total']

# =============================================================================
# Set up color
# =============================================================================

mut_freq['color'] = mut_freq['EFFECT'].map(mut_colors)

# =============================================================================
# Set up copy number color as well
# =============================================================================

cn_x_dict = dict(zip(cn_genes, np.arange(len(mut_genes),len(mut_genes)+len(cn_genes),1)))
cn_freq['x'] = cn_freq['GENE'].map(cn_x_dict)
cn_freq.loc[cn_freq['latitude'] == 'High-risk', 'x'] = cn_freq['x']+0.4

# =============================================================================
# CN
# =============================================================================

cn_freq['Total'] = cn_freq['latitude'].map(cn_risk_totals)
cn_freq['Frequency'] = cn_freq['Count']/cn_freq['Total']

# =============================================================================
# Keep directional
# =============================================================================

cn_freq = cn_freq[(cn_freq['Copy_num'] < 0)|(cn_freq['GENE'].isin([]))]
cn_freq['color'] = cn_freq['Copy_num'].map(cn_color)

# =============================================================================
# Plot
# =============================================================================

fig,ax = plt.subplots(figsize = (3,2.25))

# =============================================================================
# Mutations
# =============================================================================

idx = pd.MultiIndex.from_product([mut_genes, ['High-risk','Low-risk']], names = ['GENE','latitude'])
bottom = pd.DataFrame(index = idx, columns = ['b']).fillna(0)

for effect in list(mut_colors.keys()):
    tmp = mut_freq[mut_freq['EFFECT'] == effect].copy().set_index(['GENE','latitude'])
    tmp = tmp.merge(bottom, left_index = True, right_index = True)
    ax.bar(tmp['x'], tmp['Frequency'], width = 0.4, color = tmp['color'], bottom = tmp['b'])
    bottom = bottom.merge(tmp[['Frequency']], left_index = True, right_index = True, how = 'left').fillna(0)
    bottom['b'] = bottom['b']+bottom['Frequency']
    bottom = bottom.drop('Frequency', axis = 1)
    

# =============================================================================
# Copy number
# =============================================================================

cn_freq_1 = cn_freq[cn_freq['Copy_num'] == -1]
cn_freq_2 = cn_freq[cn_freq['Copy_num'] == -2]

cn_freq_2 = cn_freq_2.merge(cn_freq_1[['GENE','latitude','Frequency']].rename(columns = {'Frequency':'Freq_1'}), on = ['GENE','latitude'])

ax.bar(cn_freq_1['x'], cn_freq_1['Frequency'], color = cn_freq_1['color'], width = 0.4)
ax.bar(cn_freq_2['x'], cn_freq_2['Frequency'], bottom = cn_freq_2['Freq_1'], color = cn_freq_2['color'], width = 0.4)

# =============================================================================
# Aethetics
# =============================================================================

ax.set_xticks(np.arange(0.2,len(mut_genes)+len(cn_genes),1))
ax.set_xticklabels(mut_genes+cn_genes, rotation = 90, fontsize = 6)

ax.set_ylabel('Alteration frequency', fontsize = 6)

ax.tick_params(labelsize = 6)

plt.tight_layout()

fig.savefig('C:/Users/amurtha/Dropbox/Ghent M1 2019/Figures/Work from 2021/summary/AlterationFrequencyByRiskGroup.pdf')

# =============================================================================
# STATS FUNNNNNNNN
# =============================================================================

import scipy.stats as stats

m_stats = pd.DataFrame(index = mut_genes, columns = ['HR Altered','HR WT','LR Altered','LR WT','pval'])

for gene in mut_genes:
    hr = mut_freq[(mut_freq['GENE'] == gene)&(mut_freq['latitude'] == 'High-risk')].copy()
    m_stats.at[gene, 'HR Altered'] = hr['Count'].sum()
    m_stats.at[gene, 'HR WT'] = mut_risk_totals.get('High-risk') - hr['Count'].sum()
    
    lr = mut_freq[(mut_freq['GENE'] == gene)&(mut_freq['latitude'] == 'Low-risk')].copy()
    m_stats.at[gene, 'LR Altered'] = lr['Count'].sum()
    m_stats.at[gene, 'LR WT'] = mut_risk_totals.get('Low-risk') - lr['Count'].sum()
    
    for col in m_stats.columns:
        m_stats[col] = m_stats[col].astype(float)
    
    ct = [[m_stats.at[gene, 'HR Altered'],m_stats.at[gene, 'LR Altered']],
          [m_stats.at[gene, 'HR WT'],m_stats.at[gene, 'LR WT']]]
    
    s, p = stats.fisher_exact(ct)
    m_stats.at[gene,'pval'] = p
    
cn_stats = pd.DataFrame(index = cn_genes, columns = ['HR Altered','HR WT','LR Altered','LR WT','pval'])

for gene in cn_genes:
    hr = cn_freq[(cn_freq['GENE'] == gene)&(cn_freq['latitude'] == 'High-risk')].copy()
    cn_stats.at[gene, 'HR Altered'] = hr['Count'].sum()
    cn_stats.at[gene, 'HR WT'] = cn_risk_totals.get('High-risk') - hr['Count'].sum()
    
    lr = cn_freq[(cn_freq['GENE'] == gene)&(cn_freq['latitude'] == 'Low-risk')].copy()
    cn_stats.at[gene, 'LR Altered'] = lr['Count'].sum()
    cn_stats.at[gene, 'LR WT'] = cn_risk_totals.get('Low-risk') - lr['Count'].sum()
    
    for col in cn_stats.columns:
        cn_stats[col] = cn_stats[col].astype(float)
    
    ct = [[cn_stats.at[gene, 'HR Altered'],cn_stats.at[gene, 'LR Altered']],
          [cn_stats.at[gene, 'HR WT'],cn_stats.at[gene, 'LR WT']]]
    
    s, p = stats.fisher_exact(ct)
    cn_stats.at[gene,'pval'] = p